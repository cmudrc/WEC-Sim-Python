"""Actuator used by the published WaveStar nonlinear predictive controller."""

import numpy as np
from scipy.signal import cont2discrete


class WaveStarNmpcActuator:
    """Map a sampled torque request to axial PTO force.

    The controller's kinematic constants are retained separately from the
    published body's B-to-C geometry. The Simulink model uses these values in
    its torque-to-force block, even though its physical rod is slightly
    shorter at the neutral pose.
    """

    def __init__(self):
        self._command = [0.0, 0.0]
        self._torque = [0.0, 0.0]

    @staticmethod
    def moment_arm(stroke):
        """Return the controller's idealized metres-per-radian arm."""
        length = np.asarray(stroke, dtype=float) + .381408
        cosine = (length**2 - .412**2 - .2**2) / (-2 * .412 * .2)
        if (not np.isfinite(length).all() or np.any(length <= 0)
                or np.any(np.abs(cosine) >= 1)):
            raise ValueError("WaveStar NMPC stroke is outside the controller linkage")
        return .412 * .2 * np.sqrt(1 - cosine**2) / length

    def step(self, command_torque: float, stroke: float) -> float:
        """Advance one published 0.05 s actuator sample and return force N."""
        if not np.isfinite(command_torque):
            raise ValueError("WaveStar NMPC torque command must be finite")
        arm = self.moment_arm(stroke)
        command = float(np.clip(command_torque, -12, 12))
        torque = (
            .96294775 * command
            - 1.92342278 * self._command[0]
            + .96053421 * self._command[1]
            + 1.99777700 * self._torque[0]
            - .99783206 * self._torque[1]
        )
        self._command[:] = [command, self._command[0]]
        self._torque[:] = [torque, self._torque[0]]
        return float(np.clip(torque, -12, 12) / arm)


class WaveStarNmpcObserver:
    """Replay the published sampled position/velocity filter and Kalman observer.

    ``step`` takes the measured PTO stroke and the preceding torque command.
    It returns pitch, pitch speed, two radiation states, and estimated wave
    excitation moment. The observer uses the controller's idealized linkage,
    which is slightly different from the physical B-to-C rod geometry.
    """

    def __init__(self, dt: float = .05):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("WaveStar NMPC observer step must be positive")
        self.dt = float(dt)
        inertia = 1.04 + .4805
        continuous = np.array([
            [0., 1., 0., 0.],
            [-92.33 / inertia, -(-.1586 + 1.8) / inertia,
             -4.739 / inertia, -.5 / inertia],
            [0., 8., -13.59, -13.35],
            [0., 0., 8., 0.],
        ])
        input_matrix = np.array([[0.], [1 / inertia], [0.], [0.]])
        transition, input_discrete, _, _, _ = cont2discrete(
            (continuous, input_matrix, np.eye(4), np.zeros((4, 1))),
            self.dt, method="zoh",
        )
        self._transition = np.eye(5)
        self._transition[:4, :4] = transition
        self._transition[:4, 4] = input_discrete[:, 0]
        self._input = np.r_[input_discrete[:, 0], 0.]
        self._measurement = np.eye(5)[:2]
        self._process_covariance = np.diag([5e-6] * 4 + [.35])
        self._measurement_covariance = np.diag([6e-6] * 2)
        self._covariance = self._process_covariance.copy()
        self._state = np.zeros(5)
        self._previous_position = 0.
        self._previous_highpass = 0.
        self._previous_velocity = 0.

    @staticmethod
    def _position_from_stroke(stroke: float) -> float:
        length = float(stroke) + .381408
        cosine = (length**2 - .412**2 - .2**2) / (-2 * .412 * .2)
        if (not np.isfinite(length) or length <= 0
                or not np.isfinite(cosine) or abs(cosine) >= 1):
            raise ValueError("WaveStar NMPC stroke is outside the controller linkage")
        return float(np.arccos(cosine) - 1.170165)

    def step(self, stroke: float, previous_command: float) -> np.ndarray:
        """Advance one sample; ``previous_command`` is torque in N m."""
        if not np.isfinite(previous_command):
            raise ValueError("WaveStar NMPC previous command must be finite")
        position = self._position_from_stroke(stroke)
        frequency = 2 * np.pi * 30
        denominator = 2 + frequency * self.dt
        delayed_denominator = frequency * self.dt - 2
        highpass = (
            2 * frequency * (position - self._previous_position)
            - delayed_denominator * self._previous_highpass
        ) / denominator
        velocity = (
            frequency * self.dt * (highpass + self._previous_highpass)
            - delayed_denominator * self._previous_velocity
        ) / denominator
        self._previous_position = position
        self._previous_highpass = highpass
        self._previous_velocity = velocity

        predicted = (self._transition @ self._state
                     + self._input * previous_command)
        covariance = (self._transition @ self._covariance
                      @ self._transition.T + self._process_covariance)
        innovation_covariance = (
            self._measurement @ covariance.T @ self._measurement.T
            + self._measurement_covariance
        )
        gain = np.linalg.solve(
            innovation_covariance, self._measurement @ covariance.T,
        ).T
        measurement = np.array([position, velocity])
        self._state = (predicted + gain @ (
            measurement - self._measurement @ predicted
        ))
        self._covariance = covariance - gain @ self._measurement @ covariance
        # The pinned observer applies this calibration to its persistent
        # excitation state after every update, including the first sample.
        self._state[-1] -= .6979
        return self._state.copy()
