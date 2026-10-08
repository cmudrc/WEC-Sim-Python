# wec-sim-python

> **cmudrc fork status:** This is an active parity effort, not yet a complete
> wave energy converter simulator. Selected Sphere, RM3, OSWEC, and barge
> cases have paired MATLAB trajectory checks; general WEC-Sim dynamics and
> several published applications remain unsupported.
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
the supported device dynamics have paired MATLAB checks. The published
generalized-body-mode barge now has a coupled rigid and flexible regular-wave
runner for its floating three-DOF joint. Other GBM layouts and wave options
remain unsupported; see [PARITY.md](PARITY.md). The Python API is
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
Python builder covers `linear_subspace`, the single-body `floating_gbm`
regular-wave layout, and the paired two-body `floating_joint` layout. Mapped
coordinates use small-motion kinematics and fixed-axis PTOs. The floating
joint uses its own pitched-slider geometry and a relative-heave PTO; arbitrary
PTO attachment points are not part of that reduced layout.

For the published RM3 floating joint, configure the two bodies in HDF5 order:

```python
from wecsim import RegularCICWave, WEC

wec = WEC("RM3 floating joint")
float_body = wec.body("float", "path/to/rm3.h5",
                      inertia=(0, 21_306_090.66, 0))
spar = wec.body("spar", "path/to/rm3.h5",
                inertia=(0, 94_407_091.24, 0))
wec.floating_joint(float_body, spar, damping=1_200_000)
result = wec.run(RegularCICWave(2.5, 8), dt=0.1, end_time=400,
                 ramp_time=100, radiation_memory=60)
stroke = result.ptos["relative_heave"].stroke
```

The default uses implicit added mass. For comparisons with the pinned MATLAB
Simulink block, `added_mass_scheme="simulink_delay"` selects its numerical
feedback setting explicitly. `floating_joint` also accepts a joint surge
spring, PTO stiffness and equilibrium, hard stops, and convolution/FIR
radiation where the underlying paired solver supports them. The returned
`absorbed_power` counts the linear damper; spring energy exchange can be
computed from `-force * velocity`. The returned stroke is float heave minus
spar heave, so an initial body offset appears in the first stroke sample.

For the published generalized-body-mode barge, generate its HDF5 with the
Applications `Generalized_Body_Modes/hydroData/bemio.m`, then run:

```python
from wecsim import RegularWave, WEC

wec = WEC("GBM barge")
barge = wec.body("barge", "path/to/barge.h5",
                 inertia=(6.667e7, 2.167e9, 2.167e9))
wec.floating_gbm(barge)
result = wec.run(RegularWave(height=2, period=8),
                 dt=0.05, end_time=400, ramp_time=100)
surge = result.bodies["barge"].position[:, 0]
mode_displacement = result.flexible_modes["barge"].position
```

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
| `fixed_hinge` | `pm` or `pm_multi` | Hydrodynamic flap, optional fixed hydrodynamic base, pitch PTO | Directional PM excitation and radiation convolution; `pm_multi` sums independently phased seas |
| `fixed_hinge` | `spectrumImportFullDir` | Hydrodynamic flap, optional fixed hydrodynamic base, pitch PTO | Imported frequency-dependent directional spectrum and radiation convolution |
| `fixed_hinge` | `regular` | One hydrodynamic flap, optional fixed nonhydrodynamic base, pitch PTO | Regular-wave excitation and constant-frequency radiation |
| `fixed_morison` | `pm` | Stationary bodies without HDF5; body-local Cartesian Morison elements; no PTO | Directional irregular-wave velocity, acceleration, and six-component Morison force |
| `linear_subspace` | `none` or zero-heading `regular` | One hydrodynamic body with axial body-local Morison elements; pure heave, or surge/heave/pitch in regular waves; no PTO | Radiation convolution or constant-frequency radiation, relative-fluid drag, fluid inertia, and Morison added mass in the acceleration solve |
| `floating_joint` | `regular` or `regularCIC` | Two equilibrium-mass bodies, pitch inertias, relative-heave PTO; optional `body_to_body` | Constant-frequency radiation, impulse-response convolution, or sampled FIR radiation |
| `floating_joint` | `elevationImport` | Two equilibrium-mass bodies, relative-heave PTO, optional joint surge spring | Imported MAT elevation and radiation convolution |
| `floating_joint` | `none` | Two equilibrium-mass bodies, named initial coordinates and speeds, relative-heave PTO | Radiation convolution for paired free decay; sampled FIR is also available |
| `linear_subspace` | `regular`, `regularCIC`, `pm`, `jonswap`, `spectrumImport`, `elevationImport`, or `none` | Any number of six-DOF hydrodynamic bodies; named motions or 6-by-N maps; optional linear or rotational PTOs; selected mean-drift coefficients for regular waves | Constant-frequency radiation or convolution |

The published OSWEC passive-yaw cases use one yaw coordinate and a torsional
PTO. Set `passive_yaw=True` on the moving body to interpolate excitation at
the wave heading relative to its current yaw angle:

```python
from wecsim import RegularWave, WEC, WorldPoint

wec = WEC("Passive-yaw OSWEC")
flap = wec.body("flap", "oswec.h5", mass=12_700,
                inertia=(1.85e6,) * 3, passive_yaw=True)
wec.body("base", "oswec.h5", mass=999, inertia=(999,) * 3)
yaw = wec.coordinate("yaw", flap.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
wec.rotational_pto("hinge", yaw, damping=120_000)
result = wec.run(RegularWave(2.5, 8, direction=10),
                 dt=0.01, end_time=600, ramp_time=100)
angle = result.bodies["flap"].position[:, 5]
absorbed_power = result.ptos["hinge"].absorbed_power
```

The rotational PTO reports angle in radians as `stroke`, angular speed in
rad/s as `velocity`, and torque in N m as `force`. The same device accepts a
Pierson–Moskowitz sea with radiation memory:

```python
from wecsim import PMWave

result = wec.run(PMWave(2.5, 8, direction=10, seed=1),
                 dt=0.01, end_time=250, ramp_time=100,
                 radiation_memory=40)
```

`PMWave(..., phase_file="phases.csv")` replays a saved frequency-by-direction
phase matrix; a seed creates a reproducible Python realization with a
different random sequence from MATLAB. Passive yaw supports one pure-yaw
hydrodynamic body, stationary additional bodies, and independent radiation
with full-circle BEM headings. Python interpolates the wave heading
continuously. The published MATLAB regular and irregular passive-yaw inputs
hold heading coefficients until yaw changes by 0.01° and 1°, respectively.
For the published irregular PM case, set `passive_yaw_threshold=1` in the
moving body's `wec.body(...)` call. This reproduces the source's one-degree
coefficient hold and snap to a nearby tabulated BEM heading. A zero threshold
keeps continuous interpolation. The sampled setting is available for PM waves
with one incident direction; small trajectory differences can change its
update sample and accumulate over long runs.

The published `Morison_Element/morisonElement` application has a fixed
monopile and tower without HDF5 hydrodynamic bodies. Its nonlinear drag and
fluid-inertia force can be configured with a body-local point:

```python
import numpy as np
from wecsim import PMWave, WEC

wec = WEC("fixed monopile")
monopile = wec.fixed_body("monopile", center_gravity=(0, 0, -15),
                          volume=np.pi * 10**2 * 30,
                          inertia=(1.25e9, 1.25e9, .15e9))
wec.fixed_body("tower", center_gravity=(0, 0, 25),
               mass=1_031_930, inertia=(9.66e8, 9.66e8, .132e8))
wec.morison_element(monopile, point=monopile.at(0, 0, 10),
                    drag_coefficient=(1, 1, 1),
                    added_mass_coefficient=(1, 1, 1),
                    area=(300, 300, np.pi * 10**2 / 4),
                    volume=np.pi * 10**2 * 30,
                    phase_mode="matlab_shared")
sea = PMWave(2, 5, seed=5, directions=(0, 30, 90),
             spreading=(.1, .2, .7), frequency_range=(.001, 10),
             water_depth=30)
result = wec.run(sea, dt=.01, end_time=400, ramp_time=100, rho=1025)
force_and_moment = result.body_forces["monopile"]  # time × 6, N and N m
```

The pinned MATLAB function uses the first phase column for Morison force at
every heading, while wave elevation uses each heading's own phase.
`phase_mode="matlab_shared"` selects that source behavior; the default
`"directional"` uses each heading's own phase for force. Both modes retain
WEC-Sim's per-heading drag calculation. Supply a saved
three-column phase CSV to reproduce a MATLAB realization exactly.
An axial element on a moving hydrodynamic body can also contribute drag,
fluid inertia, and added mass. For the published Sphere's surge, heave, and
pitch joint in zero-heading regular waves:

```python
from wecsim import RegularWave, WEC

wec = WEC("Sphere with axial Morison element")
sphere = wec.body("sphere", "sphere.h5",
                  inertia=(20907301, 21306090.66, 37085481.11))
for axis in ("surge", "heave", "pitch"):
    wec.coordinate(axis, sphere.move(axis))
wec.morison_element(sphere, point=sphere.at(0, 0, -2),
                    drag_coefficient=(0, 0, 1),
                    added_mass_coefficient=(0, 0, 1),
                    area=(0, 0, 100), volume=20)
result = wec.run(RegularWave(1, 8), dt=.01, end_time=40, ramp_time=10,
                 initial_coordinate={"heave": 1})
morison_wrench = result.body_forces["sphere"]  # time × 6, N and N m
```

For pure heave, define only the heave coordinate. Still-water free decay
uses `NoWave()` and `radiation_memory=15`. These moving-body layouts have
no PTO; regular waves require zero heading and a body centered at x=y=0.
The three-DOF Morison force uses a proper pitch rotation and does not copy
the pinned MATLAB source function's nonorthogonal rotation. Other moving
Morison layouts, current profiles, and the normal/tangential coefficient
mode remain unsupported.

For a heave or other `linear_subspace` device, `JONSWAPWave(2.5, 8,
seed=1, gamma=3.3)` selects a JONSWAP sea. Omitting `gamma` uses WEC-Sim's
height/period-dependent value. The published Sphere MPC case uses 2.5 m and
8 s, which infers `gamma=1`; its JONSWAP spectrum therefore equals the PM
spectrum for those inputs. `phase_file` replays a saved MATLAB realization.
For a hydrodynamic body, `PMWave` and `JONSWAPWave` also accept
`frequency_range=(0.5, 1.5)` in rad/s to narrow the BEM frequency interval.
Limits outside the HDF5 range are replaced by that range's endpoints, as in
the pinned MATLAB wave class. A `water_depth` override is currently available
only for fixed Morison bodies; hydrodynamic bodies use their HDF5 depth.
The published Sphere MPC controller has a focused Python runner:

```python
from wecsim import run_sphere_mpc

result = run_sphere_mpc(
    "sphere.h5", "coeff.mat", seed=1,
    max_force=2e6, max_force_rate=1.5e6,
)
heave = result.position
pto_force = result.pto_force
```

This reproduces the published single-body heave layout, fourth-order
radiation fit, 0.5 s optimizer updates, and 0.5 s command transition. The
force, force-rate, heave, speed, horizon, and penalty settings are adjustable.
An integer seed produces a reproducible Python sea; pass a 500-by-1 `phase`
array to replay MATLAB's sea. The paired gate covers the published defaults,
not arbitrary changes to the controller or a general multi-body MPC.

For the published Sphere `Mean_Drift` application, select the control-surface
coefficient in the generated HDF5 file and use convolution radiation:

```python
from wecsim import RegularCICWave, WEC

wec = WEC("Sphere mean drift")
sphere = wec.body("sphere", "Mean_Drift/hydroData/sphere.h5",
                  inertia=(837.75804096,) * 3,
                  mean_drift="control_surface")
for axis in ("surge", "heave", "pitch"):
    wec.coordinate(axis, sphere.move(axis))
result = wec.run(RegularCICWave(0.1, 2), dt=0.01, end_time=100,
                 ramp_time=20, radiation_memory=10)
drift_force = dict(result.raw.extra_outputs)["body1_mean_drift_force"]
excitation_force = dict(result.raw.extra_outputs)["body1_excitation_force"]
```

`mean_drift="momentum_conservation"` selects that HDF5 coefficient when
present. An absent selected dataset raises an error. The regular-wave drift
force is the selected coefficient times wave amplitude squared and the wave
force ramp. The published unmoored linear Sphere run travels more than 16 m
in surge over 100 s; the paired result reproduces that source model but does
not establish physical accuracy at such a large displacement.

The runner rejects unsupported layouts and settings. For PM or JONSWAP, supply
`wave.height`, `wave.period`, optional `directions` and `spreading`, and either
an integer `seed` or a `phase_file` CSV to replay a MATLAB realization.
JONSWAP also accepts an optional positive `gamma`.
For `linear_subspace`, `wave.frequency_range` optionally narrows the HDF5
frequency interval in rad/s. A `wave.water_depth` override is supported only
for fixed Morison bodies.
`simulation` accepts `dt`, `end_time`, optional `ramp_time`, `rho`, `g`, and
`radiation_memory` for radiation-memory cases. For an RM3 `regularCIC` floating
joint, set `simulation.radiation_method` to `"fir"` to use the published
discrete FIR calculation; `"convolution"` remains the default. The same
setting is available for no-wave free decay, without a paired FIR baseline. Both use
`radiation_memory` (60 s by default). The general dynamics module
assembles the supported body and PTO forces; it does not parse Simscape models.
For a `linear_subspace` device, the Python API also accepts a WEC-Sim
three-column imported spectrum MAT file, including its saved phases:

```python
from wecsim import ImportedSpectrumWave

result = wec.run(ImportedSpectrumWave("spectrumData1.mat"),
                 dt=0.1, end_time=400, ramp_time=100,
                 radiation_memory=60, base_dir="path/to/inputs")
```

The file path resolves from `base_dir`. This selects an incident sea for the
configured Python WEC; the specialized RM3 floating-joint MCR runner remains
the paired solver for the published four-coordinate RM3 sea-state motion.
For a sampled time/elevation MAT record, use the same builder with
`ImportedElevationWave`:

```python
from wecsim import ImportedElevationWave

result = wec.run(ImportedElevationWave("etaData.mat"),
                 dt=0.01, end_time=40, ramp_time=10,
                 radiation_memory=15, base_dir="path/to/inputs")
```

The named MAT variable defaults to `etaData` and must contain increasing time
and elevation columns. `reapply_force_ramp=True` explicitly reproduces the
second force ramp used by the pinned MATLAB body block in paired comparisons.
The inherited `WaveClass` now uses the pinned MATLAB PM and JONSWAP spectrum
definitions, including height-dependent PM energy and JONSWAP's inferred
`gamma` when it is unspecified. Its seeded phases use a local NumPy generator;
MATLAB's Threefry substreams produce different realizations for the same
integer seed. To replay a MATLAB realization, assign its frequency-by-direction
phase matrix before `waveSetup`:

```python
import numpy as np
from wecsim.waveClass import WaveClass

wave = WaveClass("irregular")
wave.T, wave.H = 8, 2.5
wave.spectrumType = "PM"
wave.freqDisc, wave.numFreq = "Traditional", 64
wave.waveDir, wave.waveSpread = [0, 30, 90], [0.1, 0.2, 0.7]
wave.phaseData = np.loadtxt("pm_phase.csv", delimiter=",")
wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)
```

The historical BS fixture remains a compatibility check for the
original Python port; current MATLAB WEC-Sim rejects BS inputs.
For this RM3 convolution layout, the solver follows the published pitched
slider geometry and uses an implicit effective added mass by default. To
reproduce the pinned MATLAB/Simulink numerical trajectory, set
`simulation.added_mass_scheme` to `"simulink_delay"` (or pass that keyword to
`solve_rm3_regular` and the RM3 MCR runners). This explicitly selects the
source model's mass split and 1e-7 s acceleration delay; the delay is a
numerical setting, not a WEC property. The published Cases 5 and 6 use a
suspect fitted state-space radiation model and remain unsupported. The paired
Cases 3 and 4 tests check the ordinary implicit solver separately from this
opt-in numerical comparison; neither path uses the Cases 5 and 6 fit.
For the `fixed_hinge` layout, an optional second body can be fixed. The
regular-wave case accepts a nonhydrodynamic base declared with `nonhydro:
true`, `fixed: true`, and a three-component `center_gravity`; PM cases can
instead use a hydrodynamic base with `fixed: true`. Its stationary motion
appears in the response. With this base, `constraint.location` is its ground
attachment and `pto.location` is the flap's hinge attachment. The base's
constraint reaction forces are not yet calculated.
The `regularCIC` floating-joint path has paired MATLAB checks for RM3
body-to-body Cases 3 and 4 and all eight physical settings in the published
RM3 Multiple Condition Runs Option 1 sweep. The published RM3 radiation
options case also checks constant, convolution, and FIR dynamics over 500 s.
MATLAB's fitted radiation state-space Cases 5 and 6 remain outside that
validated path.
The programmatic RM3 solver also accepts `mooring_surge_stiffness` for a
linear spring at its floating joint and `excitation_force` for sampled imported
waves. The pinned `MooringMatrix` case pairs these forces and the full 400 s
body, PTO, and mooring trajectories with MATLAB. Its source body block applies
the wave ramp again after the imported-elevation convolution; the paired test
reproduces that source-specific step explicitly. The same setup can be passed
to the public Python case runner without precomputing a force array:

```python
from wecsim import run_case

case = {
    "simulation": {"dt": 0.01, "end_time": 400, "ramp_time": 40,
                   "radiation_memory": 60},
    "wave": {"type": "elevationImport", "file": "Mooring/MooringMatrix/etaData.mat",
             "variable": "etaData", "reapply_force_ramp": True},
    "bodies": [
        {"hydro_file": "_Common_Input_Files/RM3/hydroData/rm3.h5",
         "hydro_body": 1, "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
        {"hydro_file": "_Common_Input_Files/RM3/hydroData/rm3.h5",
         "hydro_body": 2, "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
    ],
    "constraint": {"kind": "floating_joint", "location": [0, 0, 0],
                   "initial_coordinate": {"spar_heave": -0.21}},
    "pto": {"kind": "relative_heave", "damping": 1_200_000},
    "mooring": {"kind": "joint_surge_spring", "stiffness": 100_000},
}
response = run_case(case, base_dir="path/to/WEC-Sim_Applications")
```

`reapply_force_ramp=True` reproduces the pinned MATLAB body block's second
force ramp; its default is `False`. The current mooring setting is a surge
spring at this joint. Arbitrary mooring matrices and attachment locations
remain unsupported.

The RM3 floating-joint solver accepts optional PTO hard stops in Python:

```python
from wecsim import LinearHardStops, solve_rm3_regular

stops = LinearHardStops(
    lower_bound=-0.6, upper_bound=0.6,
    lower_stiffness=100_000_000, upper_stiffness=100_000_000,
)
response = solve_rm3_regular(
    "rm3.h5", pto_hard_stops=stops, dt=0.025, end_time=400,
)
print(response.pto_stroke, response.pto_stop_force)
```

The case runner also accepts these names under `pto.hard_stops` for a
`floating_joint` with regular waves. Hard stops select adaptive integration,
constant-frequency radiation, and implicit added mass. Ordinary RM3 cases
retain their existing solver. The published MATLAB End_Stops run is
time-step sensitive after contact. Paired motion, force, and energy checks
pass through 400 s against the 0.0125 s MATLAB run; a 0.025 s run gates source
time-step convergence. See
`PARITY.md` for the limits and the published 0.1 s discrepancy.

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
Force limits, general controller networks, and hydraulic PTO models are not
implemented.

### Reactive direct-drive PTO

The published `Controls/ReactiveWithPTO` Sphere case couples PI control to a
simple direct-drive generator. Its winding inductance and resistance set the
generator torque response time; gear ratio reflects drivetrain inertia and
friction into the heave equation. Configure those properties on a PTO:

```python
from wecsim import RegularWave, SimpleDirectDrive, WEC, WorldPoint

wec = WEC("reactive sphere with direct drive")
sphere = wec.body("sphere", "sphere.h5", mass="equilibrium")
wec.coordinate("heave", sphere.move("heave"))
wec.pto(
    "PTO1", WorldPoint(0, 0, -2), sphere.at(0, 0, 0), axis=(0, 0, 1),
    direct_drive=SimpleDirectDrive(
        kp=380_000, ki=-152_000, torque_constant=7.186,
        gear_ratio=100, drivetrain_inertia=2, drivetrain_friction=1,
        winding_resistance=.483, winding_inductance=5.223e-3,
    ),
)
result = wec.run(RegularWave(height=2.5, period=9.6664),
                 dt=.01, end_time=200, ramp_time=50)
drive = result.ptos["PTO1"].direct_drive
# drive.current, drive.voltage, drive.resistance_loss,
# drive.electrical_power, drive.shaft_torque, drive.mechanical_power
```

`electrical_power` follows the pinned Simulink block's signed `voltage *
current + resistance_loss` output. `result.ptos["PTO1"].absorbed_power` is
positive when mechanical power enters the PTO. Current validation covers one
pure-heave body, zero-heading regular waves, and one vertical direct-drive
connection. The model does not impose voltage, current, or stroke limits.

### RM3 direct linear generator

The published `PTO-Sim/RM3/RM3_DD_PTO` case uses a three-phase linear
generator between the float and spar. Its flux states, electrical angle,
phase currents and voltages, friction, and generated power are available
through a focused two-heave runner:

```python
from wecsim import DirectLinearGenerator, run_rm3_direct_linear_generator

generator = DirectLinearGenerator(
    stator_resistance=4.58, friction=-100, pole_pitch=.072,
    magnet_flux=8, inductance=.285, load_resistance=-117.6471,
)
result = run_rm3_direct_linear_generator("rm3.h5", generator=generator)
print(result.absorbed_power.mean(), result.electrical_power.mean())
```

The negative load resistance and friction follow the source PTO-Sim block's
sign convention. `absorbed_power` and `electrical_power` are positive when
the device absorbs mechanical power and delivers power to the load. This
runner represents the published vertical two-body layout.

The same generator can be configured on body-local PTO endpoints in a
`WEC` with named motions:

```python
from wecsim import DirectLinearGenerator, RegularWave, WEC

generator = DirectLinearGenerator(
    stator_resistance=4.58, friction=-100, pole_pitch=.072,
    magnet_flux=8, inductance=.285, load_resistance=-117.6471,
)
wec = WEC("RM3 with linear generator")
float_body = wec.body("float", "rm3.h5")
spar = wec.body("spar", "rm3.h5")
wec.coordinate("float_heave", float_body.move("heave"))
wec.coordinate("spar_heave", spar.move("heave"))
wec.pto("PTO1", spar.at(0, 0, 0), float_body.at(0, 0, 0),
        axis=(0, 0, 1), linear_generator=generator)
result = wec.run(RegularWave(height=2.5, period=8),
                 dt=.0005, end_time=400, ramp_time=100)
phase_current = result.ptos["PTO1"].linear_generator.phase_current
```

This public configuration is currently limited to constant-radiation
regular waves and cannot be combined with sampled PTO controls. The
attachment points use the existing small-motion, fixed-axis PTO geometry.

### Instantaneous nonlinear hydrodynamics for heave

The Python builder supports the published `Nonlinear_Hydro/ode4/Regular` and
`ode4/RegularCIC` ellipsoid cases through a heave-only mesh mode. Supply the
BEM HDF5 file and an STL
whose triangle coordinates are relative to the body's center of gravity:

```python
import numpy as np
from wecsim import RegularWave, WEC, WorldPoint

wec = WEC("ellipsoid")
body = wec.body(
    "ellipsoid", "hydroData/ellipsoid.h5",
    geometry_file="geometry/elipsoid.stl",
    nonlinear_hydro="instantaneous", mass="equilibrium",
    inertia=(1_375_264, 1_375_264, 1_341_721),
    drag_coefficient=1, drag_area=np.pi * 5**2,
)
wec.coordinate("heave", body.move("heave"))
wec.pto("PTO1", WorldPoint(0, 0, -12.5), body.at(0, 0, 0),
        damping=1_200_000)
result = wec.run(RegularWave(height=4, period=6), dt=0.05,
                 end_time=150, ramp_time=50, rho=1025)
```

For convolution radiation, use `RegularCICWave(height=4, period=6)` and
`radiation_memory=60` in `wec.run`.

The mode integrates mesh buoyancy, instantaneous free-surface
Froude–Krylov correction, heave quadratic drag, BEM diffraction/radiation,
and the configured PTO. The STL determines equilibrium mass when
`mass="equilibrium"`; this can differ from the HDF5 displaced volume.
Current validation covers one pure-heave body with its center of gravity at
horizontal origin in zero-direction regular waves, with either constant or
convolution radiation.
The published ode45 variants are tracked separately: MATLAB applies mesh
buoyancy from the preceding 0.05 s sample while the Python mode evaluates it
at the current state. Their motion comparisons are recorded as solver
diagnostics in [PARITY.md](PARITY.md), not as ode45 numerical parity.
Other motions and sea states raise an error until their mesh force and dynamics
checks are paired with MATLAB.

### Variable draft and mass in heave

The published `Variable_Hydro/Variable_Mass` sphere switches among nine BEM
datasets, masses, and draft equilibria every 100 seconds. Define those states
in Python after generating `draft1.h5` through `draft9.h5` with its `bemio.m`:

```python
import numpy as np
from wecsim import HydroState, RegularWave, WEC, WorldPoint

rho = 1025
drafts = range(1, 10)
states = [
    HydroState(
        f"hydroData/draft{draft}.h5",
        mass=rho * np.pi / 3 * draft**2 * (15 - draft),
    )
    for draft in drafts
]
wec = WEC("variable mass sphere")
sphere = wec.variable_body("sphere", states, switch_times=range(100, 900, 100))
wec.coordinate("heave", sphere.move("heave"))
wec.pto("PTO1", WorldPoint(0, 0, 2), sphere.at(0, 0, 0),
        axis=(0, 0, 1), damping=200_000)
result = wec.run(RegularWave(height=1, period=8), dt=0.01,
                 end_time=900, ramp_time=0, rho=rho)
```

Each HDF5 file supplies the corresponding equilibrium center, displaced
volume, added mass, radiation damping, restoring stiffness, and excitation.
The mass and PTO settings come from the Python configuration. Validation covers
one pure-heave body, a regular zero-heading wave, and one passive vertical PTO.
The source model's large late-time excursions are outside the small-motion
range in which linear BEM coefficients can be trusted.

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

For an imported full-directional spectrum, configure the Python case with
`wave.type = "spectrumImportFullDir"`, a MAT `wave.file`, and either a
`wave.phase_file` CSV or an integer `wave.seed`. The default includes each
heading-bin width in both elevation and force, preserving the spectrum's
integrated energy. The published MATLAB `Full_Directional_Waves` force block
omits that width even though its elevation includes it. To reproduce that
specific source trajectory, set `wave.force_quadrature` to
`"matlab_omitted"` and `wave.excitation_interpolation` to
`"spline_frequency"`; `simulation.added_mass_scheme = "simulink_delay"`
is a separate opt-in Simulink numerical comparison. The ordinary defaults
remain integrated forcing and implicit added mass.

The same wave calculation is available directly in Python:

```python
from wecsim import (imported_full_directional_components,
                    synthesize_full_directional_response)

components = imported_full_directional_components(
    "path/to/oswec.h5", "path/to/fullDirSpectrum.mat", seed=7,
)
incident = synthesize_full_directional_response(
    "path/to/oswec.h5", components,
    dt=0.05, end_time=400, ramp_time=100,
)
# incident.elevation and incident.excitation_force are NumPy arrays.
```

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
