"""Actuator used by the published WaveStar nonlinear predictive controller."""

import numpy as np


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
