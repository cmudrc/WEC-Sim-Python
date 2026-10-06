# WEC-Sim-Python

> **cmudrc fork status:** This is an active parity effort, not yet a complete
> wave energy converter simulator. Focused solvers cover Sphere free decay,
> RM3 regular-wave coupled motion, and OSWEC hinge-pitch response; the
> original general runner still has no dynamics stage.
> See [PARITY.md](PARITY.md) for verified behavior, current
> MATLAB reference revision, and the remaining work.

To run the production-code parity checks with Python 3.12:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt pytest
.venv/bin/python -m compileall -q source/objects
.venv/bin/python -m pytest -q tests/test_wave_parity.py tests/test_body_io.py tests/test_oswec_standalone.py tests/test_rm3_standalone.py tests/test_reference_cli.py
```

**WEC-Sim-Python** is Sungjun Won's Python port of
[WEC-Sim](https://github.com/WEC-Sim/WEC-Sim), the MATLAB/Simulink wave energy
converter simulator. This fork is developing and checking the Python code
against MATLAB WEC-Sim while preserving the original author's work.

## Goal of WEC-Sim-Python
**WEC-Sim-Python** aims to help researchers, start-up companies, and enthusiasts without access to MATLAB in order to use the open-source code provided by NREL and Sandia lab. Also, with growing research in the field of machine learning, **WEC-Sim-Python** could be more convenient for those who develop machine learning projects utilizing Python.

## Current status

Wave generation, RM3 and OSWEC hydrodynamic input, and tested force preprocessing
have focused checks. The main runner completes preprocessing but its general
device dynamics solver has not been implemented. The published Sphere free-decay
cases can be calculated with the focused solver:

```python
from source.objects.linearHeave import solve_heave_free_decay

response = solve_heave_free_decay("path/to/sphere.h5", initial_displacement=1.0)
# response.time, response.position, response.velocity, response.force_total
```

This solver assumes one heave-only body, zero incident waves, and no PTO,
mooring, or nonlinear force. Generate `sphere.h5` with the published
WEC-Sim_Applications Sphere `bemio.m`, or download the HDF5 artifact from the
[MATLAB reference-model run](https://github.com/cmudrc/WEC-Sim-Python/actions/runs/37472556464).
The RM3 regular-wave heave subsystem can be calculated with the same module:

```python
from source.objects.linearHeave import solve_two_body_regular_heave

response = solve_two_body_regular_heave(
    "path/to/rm3.h5", wave_height=2.5, wave_period=8.0,
    pto_damping=1_200_000.0,
)
# response.position and response.velocity each have two body columns;
# response.pto_force is the PTO force acting on body 1.
```

This calculation includes only vertical translation, a linear relative-motion
PTO, and frequency-dependent hydrodynamic coefficients at the incident wave
frequency. The coupled RM3 reference model also predicts both body surge
motions and their shared pitch:

```python
from source.objects.rm3Regular import solve_rm3_regular

response = solve_rm3_regular("path/to/rm3.h5")
# response.body_position and response.body_velocity have shape
# (time_steps, 2 bodies, 6 DOFs); response.pto_force is the heave PTO force.
```

This model uses the published RM3 inertias, regular wave, floating-joint
geometry, and linear PTO. It covers the canonical case's active degrees of
freedom; it does not supply general Simscape joint dynamics.
The OSWEC reference can generate a seeded Python wave realization and pass
its six-component excitation history to the hinge-pitch solver:

```python
from source.objects.hingePitch import solve_hinged_pitch_from_excitation
from source.objects.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)

components = pm_equal_energy_components(
    "path/to/oswec.h5", significant_height=2.5, peak_period=8,
    directions=[0, 30, 90], spreading=[0.1, 0.2, 0.7], seed=7,
)
wave = synthesize_irregular_response(
    "path/to/oswec.h5", components, dt=0.1, end_time=400, ramp_time=100,
)
response = solve_hinged_pitch_from_excitation(
    "path/to/oswec.h5", wave.excitation_force, hinge_z=-8.9,
    body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
)
# response.angle is in radians; wave.elevation is the incident elevation.
```

It models pitch about a fixed hinge and the radiation-memory force. PM
equal-energy bins, directional excitation, and pitch are checked against
current MATLAB WEC-Sim using the same saved random phase matrix. A Python
integer seed creates a reproducible Python realization, with a different
random sequence from MATLAB.

To run a supported case without writing Python code:

```sh
python -m source.objects.referenceRunner rm3 --h5 path/to/rm3.h5 --output results/rm3.csv
python -m source.objects.referenceRunner rm3 --h5 path/to/rm3.h5 --output results/rm3-b2b.csv --b2b
python -m source.objects.referenceRunner oswec --h5 path/to/oswec.h5 --output results/oswec.csv --seed 7
python -m source.objects.referenceRunner sphere --h5 path/to/sphere.h5 --output results/sphere.csv --initial-displacement 1
```

Each command writes a numeric CSV and an adjacent JSON file with the model
settings, HDF5 SHA-256 hash, NumPy version, and Git revision/dirty state.
For RM3, `--b2b` includes cross-body hydrodynamic coupling as in the
published B2B Case 2; the default matches B2B Case 1.
The original README projected completion in
August 2022; that date is no longer applicable. See [PARITY.md](PARITY.md)
for the tested scope and next reference case.

## Design direction

The Python implementation will solve device dynamics without Simulink. Its
preprocessing, equations, and numerical outputs must be compared to explicit
MATLAB WEC-Sim revisions and reference cases before a feature is called
supported. ParaView file output and BEMIO conversion are future work.
