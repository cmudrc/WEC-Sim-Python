"""Directional irregular-wave excitation from WEC-Sim hydrodynamic data.

The formulas follow the current MATLAB WEC-Sim ``irregWaveSpectrum``,
``waveElevIrreg``, and ``irregExcF`` path for a PM spectrum. Callers may supply
the phase matrix to replay a MATLAB realization or generate a reproducible
Python realization with an integer seed.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from .bodyClass import BodyClass


@dataclass(frozen=True)
class IrregularComponents:
    omega: np.ndarray
    spectral_amplitude: np.ndarray
    d_omega: np.ndarray
    directions: np.ndarray
    spreading: np.ndarray
    phase: np.ndarray


@dataclass(frozen=True)
class IrregularResponse:
    time: np.ndarray
    elevation: np.ndarray
    excitation_force: np.ndarray


def imported_spectrum_components(
    h5_file: str | Path,
    spectrum_file: str | Path,
) -> IrregularComponents:
    """Read a three-column WEC-Sim spectrumImport MAT file with saved phases.

    The columns are frequency in Hz, spectral density in m²/Hz, and phase in
    radians. The current MATLAB wave class keeps only BEM-range frequencies,
    computes midpoint bin widths in rad/s, and converts density to m²/(rad/s).
    Its phase column makes the realized sea state deterministic.
    """
    source = loadmat(spectrum_file)
    if "spectrumData" not in source:
        raise ValueError("imported spectrum MAT file needs spectrumData")
    values = np.asarray(source["spectrumData"], dtype=float)
    if (values.ndim != 2 or values.shape[1] != 3
            or not np.isfinite(values).all()):
        raise ValueError("spectrumData must have three finite columns")

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    bem_omega = np.asarray(body.hydroData["simulation_parameters"]["w"]).ravel()
    if len(bem_omega) < 2 or not np.isfinite(bem_omega).all():
        raise ValueError("hydrodynamic frequency range is invalid")
    frequency = values[:, 0]
    keep = ((frequency >= bem_omega.min() / (2 * np.pi))
            & (frequency <= bem_omega.max() / (2 * np.pi)))
    selected = values[keep]
    if (len(selected) < 2 or np.any(np.diff(selected[:, 0]) <= 0)
            or np.any(selected[:, 1] < 0)):
        raise ValueError("imported BEM-range spectrum needs increasing frequencies and nonnegative density")
    omega = selected[:, 0] * (2 * np.pi)
    width = np.empty(len(omega))
    width[0] = omega[1] - omega[0]
    width[-1] = omega[-1] - omega[-2]
    width[1:-1] = (omega[2:] - omega[:-2]) / 2
    return IrregularComponents(
        omega=omega,
        spectral_amplitude=selected[:, 1] / np.pi,
        d_omega=width,
        directions=np.array([0.0]),
        spreading=np.array([1.0]),
        phase=selected[:, 2, None],
    )


def pm_equal_energy_components(
    h5_file: str | Path,
    *,
    significant_height: float,
    peak_period: float,
    directions: np.ndarray,
    spreading: np.ndarray,
    count: int = 500,
    seed: int | None = None,
    phase: np.ndarray | None = None,
) -> IrregularComponents:
    """Build current WEC-Sim PM equal-energy bins from the BEM frequency range.

    ``phase`` overrides random generation for replay. With ``seed``, NumPy's
    generator makes a reproducible *Python* realization; its random sequence
    is not MATLAB Threefry's. Exactly matching a MATLAB run requires its
    saved phase matrix.
    """
    if (not np.isfinite([significant_height, peak_period]).all()
            or significant_height <= 0 or peak_period <= 0):
        raise ValueError("significant_height and peak_period must be positive and finite")
    if not isinstance(count, int) or count < 2:
        raise ValueError("count must be an integer of at least two")
    direction, spread = _directions_and_spread(directions, spreading)
    if phase is not None and seed is not None:
        raise ValueError("supply either phase or seed")

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    bem_omega = np.asarray(body.hydroData["simulation_parameters"]["w"]).ravel()
    if len(bem_omega) < 2 or not np.isfinite(bem_omega).all():
        raise ValueError("hydrodynamic frequency range is invalid")
    # Current MATLAB EqualEnergy uses 500,000 equal-width intervals before
    # locating the closest cumulative-energy boundary for each bin.
    dense_omega = np.linspace(bem_omega.min(), bem_omega.max(), 500_001)
    frequency = dense_omega / (2 * np.pi)
    b_pm = 1.25 * (1 / peak_period)**4
    a_pm = b_pm * (significant_height / 2)**2
    spectrum_hz = a_pm * frequency**-5 * np.exp(-b_pm * frequency**-4)
    df = frequency[1] - frequency[0]
    integrated = np.empty(len(frequency))
    integrated[0] = 0.0
    integrated[1:] = np.cumsum((spectrum_hz[:-1] + spectrum_hz[1:]) * df / 2)
    energy_per_bin = integrated[-1] / (count + 1)
    boundaries = np.zeros(count + 2, dtype=int)
    for k in range(1, count + 2):
        candidate = int(np.searchsorted(integrated, k * energy_per_bin))
        candidate = min(max(candidate, boundaries[k - 1] + 1), len(integrated) - 1)
        earlier = candidate - 1
        if earlier > boundaries[k - 1] and (
            abs(integrated[earlier] - k * energy_per_bin)
            <= abs(integrated[candidate] - k * energy_per_bin)
        ):
            candidate = earlier
        # MATLAB's ``wn(k+1) = wn(k) + wna(k)`` moves one grid point beyond
        # the nearest cumulative-energy sample.
        boundaries[k] = min(candidate + 1, len(integrated) - 1)
    indices = boundaries[1:-1]
    omega = dense_omega[indices]
    d_omega = np.diff(np.r_[dense_omega[0], omega])
    spectral_amplitude = 2 * spectrum_hz[indices] / (2 * np.pi)
    if phase is None:
        phase_array = 2 * np.pi * np.random.default_rng(seed).random((count, len(direction)))
    else:
        phase_array = np.asarray(phase, dtype=float)
        if phase_array.shape != (count, len(direction)) or not np.isfinite(phase_array).all():
            raise ValueError("phase must have shape (count, number of directions)")
    return IrregularComponents(
        omega=omega, spectral_amplitude=spectral_amplitude,
        d_omega=d_omega, directions=direction, spreading=spread,
        phase=phase_array,
    )


def synthesize_irregular_response(
    h5_file: str | Path,
    components: IrregularComponents,
    *,
    dt: float,
    end_time: float,
    ramp_time: float,
    body_number: int = 1,
    rho: float = 1000.0,
    g: float = 9.81,
) -> IrregularResponse:
    """Return wave elevation and six-component excitation at uniform times."""
    if not np.isfinite([dt, end_time, ramp_time, rho, g]).all():
        raise ValueError("time and fluid parameters must be finite")
    if dt <= 0 or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0:
        raise ValueError("dt, rho, and g must be positive; times must be nonnegative")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    omega = np.asarray(components.omega, dtype=float).ravel()
    amplitude = np.asarray(components.spectral_amplitude, dtype=float).ravel()
    d_omega = np.asarray(components.d_omega, dtype=float).ravel()
    direction, spread = _directions_and_spread(
        components.directions, components.spreading,
    )
    phase = np.asarray(components.phase, dtype=float)
    if (len(omega) < 2 or amplitude.shape != omega.shape
            or d_omega.shape != omega.shape
            or phase.shape != (len(omega), len(direction))
            or not np.isfinite(omega).all() or not np.isfinite(amplitude).all()
            or not np.isfinite(d_omega).all() or not np.isfinite(phase).all()
            or np.any(np.diff(omega) <= 0) or np.any(amplitude < 0)
            or np.any(d_omega <= 0)):
        raise ValueError("wave components have inconsistent or invalid values")

    body = BodyClass(str(h5_file))
    body.bodyNumber = body_number
    body.readH5file()
    body.irrExcitation(omega, len(omega), direction, rho, g)
    real = np.transpose(body.hydroForce["fExt"]["re"], (1, 0, 2))
    imaginary = np.transpose(body.hydroForce["fExt"]["im"], (1, 0, 2))
    mean_drift = np.transpose(body.hydroForce["fExt"]["md"], (1, 0, 2))
    wave_energy = amplitude[:, None] * d_omega[:, None] * spread[None, :]
    height = np.sqrt(wave_energy)
    drift_force = np.einsum("fd,fdc->c", wave_energy, mean_drift)

    time = np.arange(steps + 1) * dt
    elevation = np.zeros(len(time))
    excitation = np.zeros((len(time), 6))
    for start in range(0, len(time), 128):
        stop = min(start + 128, len(time))
        t = time[start:stop]
        phase_angle = t[:, None, None] * omega[None, :, None] + phase[None, :, :]
        cosine = np.cos(phase_angle) * height[None, :, :]
        sine = np.sin(phase_angle) * height[None, :, :]
        elevation[start:stop] = cosine.sum(axis=(1, 2))
        excitation[start:stop] = (
            drift_force + np.einsum("tfd,fdc->tc", cosine, real)
            - np.einsum("tfd,fdc->tc", sine, imaginary)
        )
    ramp = np.ones(len(time))
    if ramp_time > 0:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    return IrregularResponse(
        time=time, elevation=ramp * elevation,
        excitation_force=ramp[:, None] * excitation,
    )


def _directions_and_spread(directions, spreading):
    direction = np.asarray(directions, dtype=float).ravel()
    spread = np.asarray(spreading, dtype=float).ravel()
    if (len(direction) == 0 or direction.shape != spread.shape
            or not np.isfinite(direction).all() or not np.isfinite(spread).all()
            or np.any(spread < 0) or not np.isclose(spread.sum(), 1.0, atol=1e-10)):
        raise ValueError("directions and spreading must be finite, same-sized, and sum to one")
    return direction, spread
