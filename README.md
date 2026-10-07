# wec-sim-python

> **cmudrc fork status:** This is an active parity effort, not yet a complete
> wave energy converter simulator. A case-driven dynamics runner now covers
> Sphere heave free decay, RM3 regular-wave coupled motion, and OSWEC
> hinged pitch. Other mechanical layouts and force models remain unsupported.
> See [PARITY.md](PARITY.md) for verified behavior, current
> MATLAB reference revision, and the remaining work.

Install from a clone with Python 3.12:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -c 'from wecsim import WEC; print(WEC)'
```

For a regular install without cloning, use
`python -m pip install git+https://github.com/cmudrc/wec-sim-python.git`.
Both methods install the `wecsim` import and the `wecsim` and
`wecsim-reference` commands.

For development, run the production-code parity checks with:

```sh
.venv/bin/python -m pip install pytest
.venv/bin/python -m compileall -q wecsim
.venv/bin/python -m pytest -q tests/test_wave_parity.py tests/test_body_io.py tests/test_oswec_standalone.py tests/test_rm3_standalone.py tests/test_reference_cli.py tests/test_general_dynamics.py tests/test_case_dynamics.py tests/test_pto_connections.py tests/test_python_api.py
```

**WEC-Sim-Python** is Sungjun Won's Python port of
[WEC-Sim](https://github.com/WEC-Sim/WEC-Sim), the MATLAB/Simulink wave energy
converter simulator. This fork is developing and checking the Python code
against MATLAB WEC-Sim while preserving the original author's work.

The installable code is in `wecsim/`, and runnable cases and bundled RM3
inputs are in `examples/`. Sungjun Won's early sketches remain available in
Git history.

## Goal of WEC-Sim-Python
**WEC-Sim-Python** aims to help researchers, start-up companies, and enthusiasts without access to MATLAB in order to use the open-source code provided by NREL and Sandia lab. Also, with growing research in the field of machine learning, **WEC-Sim-Python** could be more convenient for those who develop machine learning projects utilizing Python.

## Current status

Wave generation, RM3 and OSWEC hydrodynamic input, force preprocessing, and
the supported device dynamics have paired MATLAB checks. The Python API is
the primary way to configure a linearized device. It constructs bodies,
named motions, attachment points, PTOs, and waves as Python objects, then
returns NumPy arrays directly. Import it as `wecsim` after installation:

```python
from wecsim import NoWave, WEC, WorldPoint

wec = WEC("Heaving float")
float_body = wec.body("float", "path/to/hydro.h5")
wec.coordinate("heave", float_body.move("heave"))
wec.pto("main", float_body.at(0, 0, 0), WorldPoint(0, 0, 10),
        damping=1_200_000)
result = wec.run(NoWave(), dt=0.1, end_time=20,
                 radiation_memory=15,
                 initial_coordinate={"heave": 1})
heave = result.bodies["float"].position[:, 2]
pto_force = result.ptos["main"].force
```

The [full Python example](examples/configurable_rm3_pto.py) defines a
two-body device with a shared pitch pivot and body-local PTO points. Run it
with `python -m examples.configurable_rm3_pto`. Relative HDF5 paths in
the Python API resolve from the current directory unless `base_dir` is passed
to `wec.run`. Body order must match the HDF5 hydrodynamic body order. The
Python builder currently covers the `linear_subspace` layout with regular or
no waves; its named coordinates use small-motion kinematics and fixed-axis
PTOs. The other validated reference layouts remain available through the
case runner and focused solver functions.
For a regular-wave declutching PTO, pass
`control=DeclutchingControl(gain=232_020, declutch_time=0.8)` to `wec.pto`
instead of constant damping. Import `DeclutchingControl` from `wecsim`.
The controller disengages at a velocity sign change and reengages after the
configured interval; `minimum_on_time` defaults to 0.2 s. The published
Sphere case is paired against MATLAB at a 0.01 s step. Other geometries and
wave settings are not yet paired.
For the published Sphere latching case, pass
`control=LatchingControl(gain=49_181, latch_damping=37_308_296, latch_time=2.4)`.
The controller applies the larger damping for the timed interval after a
velocity reversal, then returns to its normal gain. It is a finite damping
force, not a rigid lock. `minimum_normal_time` defaults to 0.2 s. The control
settings serialize through `WEC.to_case` for reproducible saved cases.

The JSON runner remains available for saved and reproducible cases. It
accepts simulation, wave, body, constraint, and PTO settings. The included
RM3 example runs with the historical bundled HDF5:

```sh
python -m wecsim examples/rm3.json --output results/rm3.csv
```

The runner saves body position and velocity, wave elevation where applicable,
PTO force or torque, and a JSON provenance record beside the CSV. HDF5 paths
in the case file resolve relative to that file. For a comparison with current
MATLAB WEC-Sim, point the case at the HDF5 file from the pinned reference run.
The provenance record hashes the case, hydrodynamic files, and any replayed
phase CSV.

Supported combinations are:

| Constraint `kind` | Wave `type` | Bodies and PTO | Integration |
| --- | --- | --- | --- |
| `heave` | `none` | One equilibrium-mass body, initial heave displacement, no PTO | Radiation convolution |
| `fixed_hinge` | `pm` | One body, explicit mass and pitch inertia, pitch PTO | Directional PM excitation and radiation convolution |
| `fixed_hinge` | `regular` | One hydrodynamic flap, optional fixed nonhydrodynamic base, pitch PTO | Regular-wave excitation and constant-frequency radiation |
| `floating_joint` | `regular` or `regularCIC` | Two equilibrium-mass bodies, pitch inertias, relative-heave PTO; optional `body_to_body` | Constant-frequency radiation, impulse-response convolution, or sampled FIR radiation |
| `floating_joint` | `none` | Two equilibrium-mass bodies, named initial coordinates and speeds, relative-heave PTO | Radiation convolution for paired free decay; sampled FIR is also available |
| `linear_subspace` | `regular` or `none` | Any number of six-DOF hydrodynamic bodies; named motions or 6-by-N maps; optional linear PTO matrices or body-local PTO connections | Constant-frequency radiation or convolution |

The runner rejects unsupported layouts and settings. For the PM case, supply
`wave.height`, `wave.period`, optional `directions` and `spreading`, and either
an integer `seed` or a `phase_file` CSV to replay a MATLAB realization.
`simulation` accepts `dt`, `end_time`, optional `ramp_time`, `rho`, `g`, and
`radiation_memory` for radiation-memory cases. For an RM3 `regularCIC` floating
joint, set `simulation.radiation_method` to `"fir"` to use the published
discrete FIR calculation; `"convolution"` remains the default. The same
setting is available for no-wave free decay, without a paired FIR baseline. Both use
`radiation_memory` (60 s by default). The general dynamics module
assembles the supported body and PTO forces; it does not parse Simscape models.
For this RM3 convolution layout, the solver follows the published pitched
slider geometry and Simulink's delayed added-mass feedback. The tiny delay is
a numerical setting of the source model, not a WEC property.
For the regular-wave `fixed_hinge` layout, an optional second body can be
declared with `nonhydro: true`, `fixed: true`, and a three-component
`center_gravity`. Its stationary motion appears in the response. With this
base, `constraint.location` is the base's ground attachment and
`pto.location` is the flap's hinge attachment. The fixed base's constraint
reaction forces are not yet calculated.
The `regularCIC` floating-joint path has paired MATLAB checks for RM3
body-to-body Cases 3 and 4 and all eight physical settings in the published
RM3 Multiple Condition Runs Option 1 sweep. The published RM3 radiation
options case also checks constant, convolution, and FIR dynamics over 500 s.
MATLAB's fitted radiation state-space Cases 5 and 6 remain outside that
validated path.

RM3 multiple-condition runs can be configured in Python without a JSON input:

```python
from wecsim import mcr_grid, run_rm3_mcr

conditions = mcr_grid(
    heights=[1.5, 2.5], periods=[6, 8],
    damping_values=[1_200_000, 2_400_000],
)
result = run_rm3_mcr("rm3.h5", conditions)
power = result.power_matrix(damping=1_200_000)
print(power.periods, power.heights, power.absorbed_power)
first_trajectory = result.traces[0].response
```

`mcr_wave_statistics(path, damping_values)` reads the published Option 2
Excel grid; `mcr_mat_file(path)` reads the Option 3 MAT-file case table.
`run_mcr(conditions, simulate, averaging_start_time=...)` accepts a Python
callback for another configured WEC. The callback returns an `MCRTrace` with
time samples and signed absorbed PTO power; an optional response object keeps
the full trajectory available. The published RM3 example averages from
199.9 s through 400 s and reports positive absorbed power, the opposite of
MATLAB's signed PTO power column. A spring may return stored energy, so an
individual absorbed-power sample may be negative.

The paired MATLAB workflow executes the actual Option 1, 2, and 3
`wecSimMCR` drivers separately. For each option it compares all eight body
and PTO trajectories, average powers, and power matrices. The Option 1
physical conditions are also paired as scalar MATLAB runs. Phase-seed sweeps,
multiple PTOs, and other MCR postprocessing remain unverified.

For the published three imported-spectrum RM3 sea states, pass the MAT-file
table to the Python runner. Paired MATLAB checks cover the incident waves,
excitation, body and PTO trajectories, and absorbed power. Across the three
400 s runs, the largest position differences are 14.4 mm surge, 1.03 mm
heave, and 0.000638 rad pitch; mean absorbed powers differ by at most 29 W:

```python
from wecsim import run_rm3_spectrum_mcr

sea_states = run_rm3_spectrum_mcr(
    "rm3.h5", "RM3_MCROPT3_SeaState/mcrExample.mat",
)
print(sea_states.mean_absorbed_power)
first_wave = sea_states.wave_elevation[0]
first_float = sea_states.traces[0].response.body_position[:, 0, :]
```

The three `spectrumData*.mat` files must sit beside `mcrExample.mat`. Their
third column supplies the exact phase used for each wave component; no random
seed is needed. The runner uses the same 60 s convolution radiation model as
the validated RM3 cases 3 and 4.

For a no-wave floating joint, `constraint.initial_coordinate` and
`constraint.initial_speed` accept four values or dictionaries keyed by
`surge`, `float_heave`, `spar_heave`, and `pitch`. In the published RM3 PTO
extension examples, `float_heave: 5` or `spar_heave: -5` creates the same
initial +5 m relative PTO stroke while moving a different body.
For `linear_subspace`, the map's rows are surge, sway, heave, roll, pitch,
and yaw; its columns are independent generalized coordinates. A two-body
heave case, for example, maps the first body's heave to coordinate 1 and the
second body's heave to coordinate 2. `constraint.initial_coordinate` and
`initial_speed` set those coordinates, and `pto.damping_matrix` and
`stiffness_matrix` apply generalized linear forces. This layout assumes
small rotations and a constant coordinate map.
The matrices may be signed to represent active linear feedback. For example,
the published Sphere reactive PI controller uses a heave-coordinate PTO with
`stiffness_matrix: [[-573350]]` and `damping_matrix: [[49181]]`, giving
`force = -49181 * heave_speed + 573350 * heave_displacement` in the runner's
force convention. This published case has no PTO stroke limit; the paired
trajectory comparison does not validate hardware feasibility.

### Configuring body motions and PTO attachments

The [configurable two-body example](examples/configurable_rm3_pto.json)
shows a WEC definition with named bodies, independent motion coordinates, and
a PTO connection. Run it with:

```sh
python -m wecsim examples/configurable_rm3_pto.json --output results/configurable-rm3.csv
```

In `constraint.coordinates`, each named coordinate lists the body motions it
drives. Assigning the same coordinate to two bodies gives them a shared motion;
separate coordinates let them move independently. A motion names one of surge,
sway, heave, roll, pitch, or yaw and may include a `scale` (default `1`).
For a rotational motion, `pivot` can specify `{"world": [x, y, z]}` or
`{"point": [x, y, z]}` in that body's local center-of-gravity frame. The
rotation then moves the body center about that point. With no `pivot`, the
body rotates about its center of gravity. The example's shared pitch rotates
both RM3 bodies about world `[0, 0, 0]`.
`initial_coordinate` and `initial_speed` may be vectors or dictionaries keyed
by coordinate name. The original 6-by-N `coordinate_map` remains available for
advanced cases. These maps describe constant, small-motion kinematics, not an
arbitrary multibody joint solver.

`ptos` is a list of linear actuators. Each PTO has a unique `name`, `from` and
`to` endpoints, and damping and/or stiffness. A body endpoint uses
`{"body": "float", "point": [x, y, z]}`: the point is in that body's local
frame relative to its HDF5 center of gravity. A fixed world anchor uses
`{"ground": [x, y, z]}`. `axis` gives a fixed world-space force direction;
if omitted, the axis follows the line between endpoints **at the reference
pose** and then stays fixed. Stroke is projected displacement along that axis
from the reference pose. Moving an attachment point changes its moment arm
when the body rotates. The CSV includes each PTO's stroke, velocity, force,
and power absorbed by its damper, as well as named coordinate motion. This
force is positive along the axis on the `to` endpoint and opposite on `from`.
The reported absorbed power is damping times stroke velocity squared for a
constant damper and zero while a declutching PTO is disengaged. This
connection model is linearized for small rotations; it does not update the
actuator's axis as its endpoints move. Use either `ptos` or the older `pto`
matrix in one case.
In the example, each body-local PTO point has a vertical coordinate that
places it at world `z = 0` in the reference pose; its horizontal offset sets
the pitch moment arm.

### PTO tuning

The `fixed_hinge` and `floating_joint` layouts accept scalar `pto.damping`
and optional `pto.stiffness`. They also accept either
`pto.equilibrium_position` or `pto.pretension`. For hinge pitch, the position
is an angle in radians and forces are torques; for the floating joint, it is
relative heave in meters and forces are newtons. The implemented law is
`F = -damping * velocity - stiffness * (position - equilibrium_position)`.
As in MATLAB WEC-Sim, `pretension` sets the equivalent equilibrium position
to `-pretension / stiffness`; nonzero pretension or equilibrium position
requires positive stiffness. These fields default to zero, preserving the
paired reference cases. For example, add
`"stiffness": 100000, "equilibrium_position": 0.1` to the RM3 example's
`pto` object to shift its neutral relative heave by 0.1 m.

The `linear_subspace` layout uses `pto.stiffness_matrix` and
`pto.damping_matrix`, with optional `pto.equilibrium_coordinate` (N values).
Its force law is `F = -K @ (q - q_eq) - C @ q_dot`. A nonzero offset must
produce a spring force. The zero-offset cases have paired MATLAB checks;
nonzero offsets have analytical and case-level tests. Each actuator in `ptos`
can instead set `damping`, `stiffness`, and either `equilibrium_position` or
`pretension`. Connection geometry and nonzero offsets have case-level and
analytical checks. A paired MATLAB Sphere case covers nonzero stiffness,
extra damping while specifying a body-local attachment shifted 1 m in x.
Its heave-only motion cannot validate attachment-location dynamics or a
rotational moment arm.
Force limits, hard stops, controllers other than the published declutching and
reactive PI laws, and hydraulic PTO models are not implemented.

The published Sphere free-decay cases can also be calculated with the focused solver:

```python
from wecsim.linearHeave import solve_heave_free_decay

response = solve_heave_free_decay("path/to/sphere.h5", initial_displacement=1.0)
# response.time, response.position, response.velocity, response.force_total
```

This solver assumes one heave-only body, zero incident waves, and no PTO,
mooring, or nonlinear force. Generate `sphere.h5` with the published
WEC-Sim_Applications Sphere `bemio.m`, or download the HDF5 artifact from the
[MATLAB reference-model run](https://github.com/cmudrc/wec-sim-python/actions/runs/37472556464).
The RM3 regular-wave heave subsystem can be calculated with the same module:

```python
from wecsim.linearHeave import solve_two_body_regular_heave

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
from wecsim.rm3Regular import solve_rm3_regular

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
from wecsim.hingePitch import solve_hinged_pitch_from_excitation
from wecsim.irregularWave import (
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
python -m wecsim.reference rm3 --h5 path/to/rm3.h5 --output results/rm3.csv
python -m wecsim.reference rm3 --h5 path/to/rm3.h5 --output results/rm3-b2b.csv --b2b
python -m wecsim.reference oswec --h5 path/to/oswec.h5 --output results/oswec.csv --seed 7
python -m wecsim.reference sphere --h5 path/to/sphere.h5 --output results/sphere.csv --initial-displacement 1
```

Each command writes a numeric CSV and an adjacent JSON file with the model
settings, HDF5 SHA-256 hash, NumPy version, and Git revision/dirty state.
For RM3, `--b2b` includes cross-body hydrodynamic coupling as in the
published B2B Case 2; the default matches B2B Case 1.
The original README projected completion in
August 2022; that date is no longer applicable. See [PARITY.md](PARITY.md)
for the tested scope and next reference case.

## Design direction

The Python dynamics engine uses independent coordinates for each supported
constraint layout, avoiding numerical joint drift. New layouts and force
models need explicit MATLAB comparisons before they are called supported.
ParaView file output and BEMIO conversion remain future work.
