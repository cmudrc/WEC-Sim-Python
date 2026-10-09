"""Wave-surface calculations and VTP output for ParaView.

The upstream file contained unfinished MATLAB syntax and could not be imported.
Body visualization still needs a separate implementation.
"""

from pathlib import Path
from xml.etree import ElementTree

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

    def write_paraview_vtp_wave(
        self, times, directory, *, domain_size,
        num_points_x=None, num_points_y=None, water_depth=None,
        num_moordyn=0,
    ) -> tuple[Path, ...]:
        """Write a WEC-Sim-style wave grid as numbered VTP PolyData frames.

        The grid covers ``[-domain_size, domain_size]`` on both axes. This
        writes the wave and ``ground.txt`` only; body and mooring meshes are
        separate visualization outputs.
        """
        times = np.asarray(times, dtype=float)
        if (times.ndim != 1 or not len(times)
                or not np.isfinite(times).all()
                or np.any(np.diff(times) <= 0)):
            raise ValueError("wave VTP times must be finite and increasing")
        nx = self.waves.viz["numPointsX"] if num_points_x is None else num_points_x
        ny = self.waves.viz["numPointsY"] if num_points_y is None else num_points_y
        if (isinstance(nx, bool) or not isinstance(nx, (int, np.integer))
                or nx < 2 or isinstance(ny, bool)
                or not isinstance(ny, (int, np.integer)) or ny < 2):
            raise ValueError("wave VTP grid needs at least two points per axis")
        if (not np.isfinite(domain_size) or domain_size <= 0
                or isinstance(num_moordyn, bool)
                or not isinstance(num_moordyn, (int, np.integer))
                or num_moordyn < 0):
            raise ValueError("wave VTP domain and mooring count must be valid")
        depth = self.waves.waterDepth if water_depth is None else water_depth
        depth = np.asarray(depth, dtype=float)
        if depth.size != 1 or not np.isfinite(depth).all() or depth.item() <= 0:
            raise ValueError("wave VTP needs a positive water depth")

        directory = Path(directory)
        wave_dir = directory / "waves"
        wave_dir.mkdir(parents=True, exist_ok=True)
        (directory / "ground.txt").write_text(
            f"{domain_size:g}\n{depth.item():g}\n{num_moordyn}\n",
            encoding="utf-8",
        )
        x = np.linspace(-domain_size, domain_size, nx)
        y = np.linspace(-domain_size, domain_size, ny)
        X, Y = np.meshgrid(x, y)
        connectivity = np.array([
            (row * nx + col, row * nx + col + 1,
             (row + 1) * nx + col + 1, (row + 1) * nx + col)
            for row in range(ny - 1) for col in range(nx - 1)
        ])
        paths = []
        for index, time in enumerate(times, start=1):
            Z = self.waveElevationGrid(float(time), X, Y)
            if Z.shape != X.shape or not np.isfinite(Z).all():
                raise ValueError("wave VTP elevation must be finite on the grid")
            root = ElementTree.Element("VTKFile", type="PolyData", version="0.1")
            polydata = ElementTree.SubElement(root, "PolyData")
            piece = ElementTree.SubElement(
                polydata, "Piece", NumberOfPoints=str(nx * ny),
                NumberOfPolys=str((nx - 1) * (ny - 1)),
            )
            points = ElementTree.SubElement(piece, "Points")
            point_values = ElementTree.SubElement(
                points, "DataArray", type="Float32",
                NumberOfComponents="3", format="ascii",
            )
            point_values.text = "\n" + "\n".join(
                f"{xx:.5f} {yy:.5f} {zz:.5f}"
                for xx, yy, zz in zip(X.ravel(), Y.ravel(), Z.ravel())
            ) + "\n"
            polys = ElementTree.SubElement(piece, "Polys")
            cells = ElementTree.SubElement(
                polys, "DataArray", type="Int32", Name="connectivity",
                format="ascii",
            )
            cells.text = "\n" + "\n".join(
                " ".join(map(str, cell)) for cell in connectivity
            ) + "\n"
            offsets = ElementTree.SubElement(
                polys, "DataArray", type="Int32", Name="offsets",
                format="ascii",
            )
            offsets.text = " " + " ".join(
                str(4 * cell) for cell in range(1, len(connectivity) + 1)
            ) + " "
            path = wave_dir / f"waves_{index}.vtp"
            ElementTree.ElementTree(root).write(
                path, encoding="utf-8", xml_declaration=True,
            )
            paths.append(path)
        return tuple(paths)

    def write_paraview_vtp(self, *args, **kwargs):
        raise NotImplementedError("body ParaView output is not implemented")
