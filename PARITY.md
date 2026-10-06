# MATLAB WEC-Sim comparison

This fork started from THREDgroup/WEC-Sim-Python `c04ef427` (October 2021).
The current MATLAB reference reviewed for this baseline is
[WEC-Sim/WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`](https://github.com/WEC-Sim/WEC-Sim/tree/0753b2e47f2457c078751dcfe5d251d1767b80ab).
The original Python repository includes MATLAB-generated wave fixtures, but
there is no MATLAB or Octave executable on the machine used for this baseline.
The current MATLAB source was inspected; it was not executed locally. The
`MATLAB reference parity` GitHub Actions workflow checks out the pinned MATLAB
revision, runs its `waveClass` under MATLAB R2025b, and compares the generated
outputs with the production Python code. It also saves the MATLAB outputs as
a workflow artifact. Its result should be checked separately from the local
fixture tests.

| Behavior | Reference | Baseline result |
| --- | --- | --- |
| Regular-wave elevation, ramp, and three gauge positions | Original MATLAB-generated `regular_1_test` files | Production `WaveClass` agrees within `1e-12` absolute error. |
| Finite-depth wave number | Current MATLAB `calcWaveNumber.m` dispersion relation | Python residual is checked to relative tolerance `1e-13`. |
| Finite-depth regular-wave power | Current MATLAB `waveClass.m` group-velocity equation | Corrected denominator to `sinh(2kh)`; tested at 25 m depth. |
| Irregular Bretschneider equal-energy binning | Original MATLAB-era test constants | Production `WaveClass` agrees after NumPy/SciPy API updates. This is historical compatibility only: current MATLAB WEC-Sim rejects the `BS` option. |
| Wave-surface grid for no-wave, regular, and directional irregular waves | Current MATLAB `waveClass.waveElevationGrid` equations | Implemented and checked for regular and two-direction irregular cases. |
| RM3 HDF5 hydrodynamic input | Original `rm3.h5` | Both bodies load, including their names and water depth, under NumPy 2. |
| RM3 regular-wave force preprocessing | Original MATLAB-generated `body_1_test` constants and files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF; the RM3 runner now completes preprocessing. |
| End-to-end RM3 motion and PTO response | Current MATLAB RM3 example | **Not yet available.** The inherited runner has no dynamics solver; its simulation stage is commented out. |

Source comparisons: [MATLAB wave class](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/waveClass.m),
[MATLAB wave-number function](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/functions/BEMIO/calcWaveNumber.m),
and [MATLAB body class](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/bodyClass.m).

## Known differences and next reference case

Current MATLAB WEC-Sim has changed its PM/JS spectra, seeded phase generator,
object properties, and some wave inputs since the Python port was written.
The historical fixtures therefore cannot establish parity for those modes.
The next end-to-end target is the included RM3 hydrodynamic input with a
regular wave, checking body displacement and PTO force against a MATLAB run
with the same inputs and solver settings. A comparison must record both code
revisions, the HDF5 input, time step, outputs, and numerical tolerances.

The inherited `paraviewClass.py` was invalid Python mixed with unfinished
MATLAB code. The wave-surface calculation is now importable; VTP serialization
and body visualization explicitly raise `NotImplementedError` until ported.
The older object tests import duplicate copies of classes inside test folders.
The parity tests here import the production files instead.
