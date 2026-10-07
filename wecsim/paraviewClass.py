"""Wave-surface calculations for ParaView output.

The upstream file contained unfinished MATLAB syntax and could not be imported.
VTP serialization and body visualization still need a separate implementation.
"""

import numpy as np


class ParaviewClass:
    def __init__(self, waves):
        self.waves = waves

    def waveElevationGrid(self, t, X, Y):
        """Return the undisturbed wave elevation over an X/Y grid at time t.

        This follows ``waveClass.waveElevationGrid`` in MATLAB WEC-Sim for
        no-wave, regular, and frequency-resolved irregular waves. The arrays
        X and Y must have matching or broadcastable shapes.
        """
        X, Y = np.broadcast_arrays(np.asarray(X, dtype=float),
                                   np.asarray(Y, dtype=float))
        waves = self.waves
        if waves.wType in ("noWave", "noWaveCIC"):
            return np.zeros_like(X)

        if waves.wType in ("regular", "regularCIC"):
            direction = np.deg2rad(np.asarray(waves.waveDir).item())
            position = X * np.cos(direction) + Y * np.sin(direction)
            return np.asarray(waves.A).item() * np.cos(
                -np.asarray(waves.k).item() * position
                + np.asarray(waves.w).item() * t
            )

        if waves.wType in ("irregular", "spectrumImport"):
            frequencies = np.atleast_1d(waves.w)
            wavenumbers = np.atleast_1d(waves.k)
            amplitudes = np.atleast_1d(waves.A)
            intervals = np.atleast_1d(waves.dw)
            directions = np.atleast_1d(waves.waveDir)
            spreads = np.atleast_1d(waves.waveSpread)
            phases = np.asarray(waves.phase)
            if phases.shape != (len(directions), len(frequencies)):
                raise ValueError("phase must have one row per direction and one column per frequency")
            if not (len(frequencies) == len(wavenumbers) == len(amplitudes) == len(intervals)):
                raise ValueError("frequency, wavenumber, amplitude, and interval lengths must match")
            if len(directions) != len(spreads):
                raise ValueError("direction and spread lengths must match")

            elevation = np.zeros_like(X)
            for direction_index, direction in enumerate(directions):
                angle = np.deg2rad(direction)
                position = X * np.cos(angle) + Y * np.sin(angle)
                for frequency_index, frequency in enumerate(frequencies):
                    amplitude = np.sqrt(
                        amplitudes[frequency_index]
                        * intervals[frequency_index]
                        * spreads[direction_index]
                    )
                    elevation += amplitude * np.cos(
                        -wavenumbers[frequency_index] * position
                        + frequency * t
                        + phases[direction_index, frequency_index]
                    )
            return elevation

        raise NotImplementedError(
            f"wave-surface visualization is not implemented for {waves.wType!r}"
        )

    def write_paraview_vtp_wave(self, *args, **kwargs):
        raise NotImplementedError("ParaView VTP serialization is not implemented")

    def write_paraview_vtp(self, *args, **kwargs):
        raise NotImplementedError("body ParaView output is not implemented")
