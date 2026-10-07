# MATLAB WEC-Sim comparison

This fork started from THREDgroup/WEC-Sim-Python `c04ef427` (October 2021).
The current MATLAB reference reviewed for this baseline is
[WEC-Sim/WEC-Sim `0753b2e47f2457c078751dcfe5d251d1767b80ab`](https://github.com/WEC-Sim/WEC-Sim/tree/0753b2e47f2457c078751dcfe5d251d1767b80ab).
The original Python repository includes MATLAB-generated fixtures, but there
is no MATLAB or Octave executable on the local machine. The `MATLAB reference
parity` GitHub Actions workflow checks out the pinned MATLAB revision, runs
its `waveClass` under MATLAB R2025b, and compares the generated outputs with
the production Python code. The live wave comparison passed on 6 October
2026. The workflow saves the MATLAB outputs as an artifact.

| Behavior | Reference | Baseline result |
| --- | --- | --- |
| Regular-wave elevation, ramp, and three gauge positions | Original MATLAB-generated `regular_1_test` files | Production `WaveClass` agrees within `1e-12` absolute error. |
| Finite-depth wave number | Current MATLAB `calcWaveNumber.m` dispersion relation | Python residual is checked to relative tolerance `1e-13`. |
| Finite-depth regular-wave power | Current MATLAB `waveClass.m` group-velocity equation | Corrected denominator to `sinh(2kh)`; tested at 25 m depth. |
| Irregular Bretschneider equal-energy binning | Original MATLAB-era test constants | Production `WaveClass` agrees after NumPy/SciPy API updates. This is historical compatibility only: current MATLAB WEC-Sim rejects the `BS` option. |
| Wave-surface grid for no-wave, regular, and directional irregular waves | Current MATLAB `waveClass.waveElevationGrid` equations | Implemented and checked for regular and two-direction irregular cases. |
| RM3 HDF5 hydrodynamic input | Original `rm3.h5` | Both bodies load, including their names and water depth, under NumPy 2. |
| Current RM3 and OSWEC HDF5 hydrodynamic inputs | Pinned current MATLAB core examples | Both bodies in each model load and complete regularCIC or directional irregular preprocessing with finite restoring, added-mass, excitation, and radiation arrays. This is input compatibility, not motion parity. |
| RM3 regular-wave force preprocessing | Original MATLAB-generated `body_1_test` constants and files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF; the RM3 runner now completes preprocessing. |
| RM3 irregular-wave force preprocessing | Original MATLAB-generated `body_2_test` files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF. |
| RM3 body interaction force preprocessing | Original MATLAB-generated `body_4_test` through `body_9_test` files | All six regular/regularCIC variants, including body-to-body coupling on/off and state-space radiation on/off, agree for both RM3 bodies on restoring stiffness, added mass, excitation, radiation IRF, and state-space matrices where present. |
| OSWEC directional irregular force preprocessing | Original MATLAB-generated `body_3_test` files | Production `BodyClass` agrees for restoring stiffness, added mass, radiation IRF, and three-direction excitation after replacing removed SciPy `interp2d`. |
| Sphere `noWaveCIC` heave free decay (0 m, 1 m, 1 m with Morison elements, 3 m, 5 m) | Current MATLAB WEC-Sim and MATLAB-generated Sphere HDF5 | The focused Python linear heave solver agrees over 40 s to maximum differences of 0.45 mm position, 0.64 mm/s velocity, and 363 N total force in the 5 m case. The Morison element has only x-direction coefficients in the published 1 m case, so it does not affect heave. |
| RM3 regular-wave heave, excitation, and PTO | Current MATLAB RM3 example and current RM3 HDF5 | A focused two-body linear heave solver with relative-motion PTO damping agrees locally over 400 s within 5.9 mm position and 6.1 mm/s velocity for both bodies. Heave excitation matches to less than `1e-5` N; PTO internal force differs by at most 8.2 kN over a 1.63 MN range. This reduced model does not cover surge, pitch, or full Simscape joint mechanics. |
| RM3 coupled surge, heave, pitch, and PTO | Current MATLAB RM3 example and current RM3 HDF5 | A four-coordinate two-body model uses the published pitched-slider geometry, body inertias, regular-wave forcing, and nonlinear rotation kinematics. Over 400 s, the two body surge positions differ by at most 36.8 mm, heaves by 4.18 mm, shared pitch by 0.000672 rad, and PTO force by 4.74 kN. It covers the canonical active DOFs but not general Simscape joint mechanics or other RM3 cases. |
| RM3 body-to-body Cases 1 and 2 | Pinned MATLAB Applications inputs and RM3 HDF5 generated with current BEMIO | Both published regular-wave cases run for 400 s with coupling off and on. Across the two, body positions differ by at most 36.8 mm surge, 4.57 mm heave, and 0.000948 rad pitch; PTO force differs by at most 5.38 kN over a 1.65 MN range. Turning on cross-body coupling substantially reduces the body-1 heave error against Case 2. These are reduced-model comparisons, not full Simscape mechanics. |
| RM3 body-to-body Cases 3 and 4 | Pinned MATLAB Applications `regularCIC` cases, with coupling off and on | Both the ordinary implicit-mass solver and the opt-in `simulink_delay` comparison are paired over 400 s. The ordinary solver differs by at most 67.2 mm surge, 3.64 mm heave, 0.00162 rad pitch, and 4.78 kN PTO force across both bodies and cases. The source-numerics setting reduces those maxima to 28.2 mm, 0.253 mm, 0.0000325 rad, and 0.292 kN; PTO stroke, speed, and power then differ by at most 0.235 mm, 0.244 mm/s, and 0.392 kW. Cases 5 and 6 use MATLAB's fitted radiation state-space model and are not validated by this convolution solver. |
| RM3 End_Stops force law and refined-step trajectory | Pinned MATLAB Applications `End_Stops` input; R2025b runs at published 0.1 s and diagnostic 0.025/0.0125 s steps through the published 400 s duration | Effective stroke bounds are ±0.6 m, each spring is 100 MN/m, stop damping is zero, and the 0.0001 m transition is active. The source's logged stop force matches the Python smooth-stop law within `1e-5` N even at saved samples inside the transition. The 0.025-to-0.0125 s source stroke difference is at most 2.79 mm over 400 s; the published 0.1 s run differs from the refined source by up to 286 mm after contact. The opt-in Python adaptive solver with implicit added mass agrees with the 0.0125 s MATLAB trajectory over 400 s within 0.986 mm PTO stroke, 6.77 mm/s PTO speed, 98.2 kN PTO force, and 0.139% ordinary damper energy. Paired gates cover both active body positions and velocities, PTO motion and force, and integrated energy; a separate source convergence gate checks both refined MATLAB steps. |
| RM3 radiation force options: constant, convolution, and FIR | Pinned MATLAB Applications `Radiation_Force_Options`, BEMIO-generated RM3 HDF5, 12 s regular waves | All three published 500 s settings are checked against the Python floating-joint solver. The FIR path uses sampled kernel taps with a held radiation force during each RK4 step. Across the three cases, the largest differences are 29.6 mm surge, 2.11 mm heave, 0.00140 rad pitch, 3.25 mm PTO stroke, 1.55 mm/s PTO speed, 1.86 kN PTO force, and 1.57 kW PTO power. Direct FIR radiation-force differences remain under their paired gates. MATLAB FIR differs from MATLAB convolution by up to 345 mm surge, so the methods are not interchangeable for this case. The application's state-space setting is preserved only in closed, unmerged diagnostic PR #6 because of negative low-frequency damping in its fit. |
| RM3 PTO extension float and spar free decays | Pinned MATLAB Applications `RM3_PTO_Extension` inputs and BEMIO-generated RM3 HDF5 | The no-wave floating-joint solver initializes the float at +5 m or the spar at −5 m, reproducing each published case's +5 m PTO stroke. Over 30 s, maximum active-body heave differences are 47.2 mm for float and 6.85 mm for spar; velocity differences are 86.1 mm/s and 1.87 mm/s. Both initial poses and the passive body's heave agree within numerical precision; PTO force and power are zero in Python and within `1.7e-8` in MATLAB output. |
| RM3 Multiple Condition Runs Option 1 physical sweep | Pinned MATLAB Applications RM3 input, expanded to eight scalar height × period × PTO-damping combinations | The Python `regularCIC` floating-joint solver agrees across all eight 400 s conditions within 75.1 mm surge, 0.387 mm heave, 0.000183 rad pitch, 0.413 mm PTO stroke, 0.231 mm/s PTO speed, 0.479 kN PTO force, and 0.487 kW PTO power. These scalar runs validate each physical condition; all three MCR drivers are separately checked below. |
| RM3 Multiple Condition Runs Options 1–3 orchestration | Pinned MATLAB Applications arrays, Option 2 Excel grid, Option 3 MAT file, and three actual MATLAB R2025b `wecSimMCR` runs | The Python loaders produce the same ordered eight-case table for all three published inputs. Each driver has a separate paired gate for every body's 4,001 position and velocity samples, every PTO stroke, speed, force, and power trace, the eight mean powers, and two 2×2 power matrices with period rows and height columns. Python mean absorbed powers differ from MATLAB by at most 0.155 kW or 0.059% across eight conditions in the MAT driver. The Option 2 Excel and Option 3 MAT driver outputs agree numerically exactly; Option 1 array and Option 3 MAT body positions differ by at most `1.9e-12` m. MATLAB's PTO power sign is opposite Python's positive absorbed-power convention. Phase-seed sweeps, multiple PTOs, and other postprocessing are not yet paired. |
| RM3 imported-spectrum sea-state MCR | Pinned MATLAB Applications `RM3_MCROPT3_SeaState` and three actual MATLAB MCR runs | All three imported component tables agree to `4.9e-15`, wave elevations to `1.0e-13` m, and active excitation forces and moments to `3.8e-8` N or N m. On MATLAB body velocities, the 60 s radiation convolution matches every logged force component within `5.4e-9` N or N m. The applied added-mass force reconstructs from the exported matrix and delayed acceleration within `7e-7` N or N m. The exact pitched-slider geometry and source-compatible mass feedback bring all three 400 s paired trajectories within 14.4 mm surge, 1.03 mm heave, and 0.000638 rad pitch; velocity differences are below 3.26 mm/s surge, 0.883 mm/s heave, and 0.000256 rad/s pitch. PTO stroke, speed, force, and absorbed power differ by at most 2.44 mm, 1.41 mm/s, 1.69 kN, and 1.19 kW. Mean absorbed power differs by at most 29 W. Explicit paired gates passed; PR #18 was merged. |
| Sphere passive controller | Pinned MATLAB Applications `Controls/Passive (P)` and MATLAB-generated Sphere HDF5 | The Python heave model represents the published proportional controller as an 860,870 N s/m damper. Maximum differences over 400 s are 0.063 mm position, 0.065 mm/s velocity, 56 N controller force, and 32 W controller power after accounting for MATLAB's opposite force and power signs. |
| Sphere reactive PI controller | Pinned MATLAB Applications `Controls/Reactive (PI)` and MATLAB-generated Sphere HDF5 | The published active feedback maps to signed linear PTO stiffness and damping in the existing Python heave coordinate. Against fresh MATLAB output over 200 s, maximum differences are 8.45 mm heave, 5.49 mm/s velocity, 4.85 kN controller force, and 35.9 kW controller power. MATLAB's logged force and power satisfy the published gain equations to numerical precision. The 22.74 m peak-to-peak heave is the published unconstrained linear controller's response; trajectory agreement does not establish a feasible PTO stroke or physical design. |
| Sphere declutching controller | Pinned MATLAB Applications `Controls/Declutching` and MATLAB-generated Sphere HDF5 | A configurable body-local PTO endpoint uses the published 232,020 N s/m engaged damping and switches off for 0.8 s after velocity reversal. In 40,001 paired samples over 400 s, heave differs by at most 1.99 mm and velocity by 4.89 mm/s. Both runs have 83 disengagement events, each 0.8 s long; event times differ by at most one 0.01 s sample, affecting four force samples. Where both runs are engaged, force differs by at most 1.14 kN and power by 1.06 kW; both report zero force while disengaged. Absorbed energy differs by 10.5 kJ over a 22.4 MJ MATLAB total. A full-sample force maximum is 106 kN because one-step mode differences occur at two transitions, so force parity is gated by event timing and within-mode values. |
| Sphere latching controller | Pinned MATLAB Applications `Controls/Latching` and MATLAB-generated Sphere HDF5 | The public `LatchingControl` applies the published 49,181 N s/m normal damping and 37,308,296 N s/m latch damping for 2.4 s after velocity reversal. Both force laws dissipate energy. In 40,001 paired samples over 400 s, MATLAB and Python each have 83 latch starts; start times differ by at most 0.05 s, and each completed latch lasts 2.41 s in the sampled source implementation. Maximum heave and velocity differences are 0.379 m and 0.573 m/s; the heave extrema differ by at most 0.033 m. Absorbed energy differs by 2.43 MJ against MATLAB's 299.98 MJ (0.81%). This switched high-gain case has 553 samples in different modes, so pointwise controller force is not a useful single parity metric; the paired gate checks the exact force and power laws, event timing, motion, envelope, and integrated energy. The published linear model does not establish a feasible physical stroke or rigid latch. |
| Configured Sphere spring, damper, and PTO attachment | Derived from the pinned Sphere passive case with 50,000 N/m PTO stiffness, 100,000 N s/m native PTO damping, and the PTO moved to x = 1 m | A fresh MATLAB run and the Python `WEC` builder agree within 0.051 mm position and PTO stroke, 0.054 mm/s velocity and PTO speed, 52 N combined force, and 27 W combined power. The x offset has no effect on heave-only motion, so this validates the force settings but not attachment-location dynamics. |
| OSWEC PM equal-energy waves, directional excitation, pitch, and PTO | Current MATLAB OSWEC example and current OSWEC HDF5 | Python recreates all 500 PM equal-energy bins from the HDF5 range, then synthesizes wave elevation and six-component excitation using the saved MATLAB phase matrix. Local component differences are below `6e-15`, elevation below `2e-13` m, and excitation below `1e-7` N. The fixed-hinge solver's pitch differs by at most 0.0021 rad over 400 s; PTO torque from the paired-force check differs by at most 37 N m. A Python seed produces a separate reproducible realization. This model covers pitch about a fixed hinge, not a general six-DOF device. |
| OSWEC fixed nonhydrodynamic base | Pinned MATLAB Applications `Nonhydro_Body` case and its BEMIO-generated OSWEC HDF5 | The regular-wave solver reports the stationary base and nonlinear flap motion about the PTO hinge. Against a fresh 400 s MATLAB R2025b run (4,001 samples), maximum flap position differences are 20.8 mm surge, 7.9 mm heave, and 0.00444 rad pitch; velocity differences are 15.6 mm/s surge, 6.8 mm/s heave, and 0.00328 rad/s pitch. All six excitation-force components differ by less than 1 N, the base position and velocity agree exactly, and zero PTO torque differs only by MATLAB numerical noise below `3.3e-7` N m. The base's ground-constraint reaction forces are not calculated. |
| Case-driven dynamics runner | Current MATLAB RM3 and OSWEC examples plus Sphere and RM3 body-to-body Applications cases | One generalized-coordinate engine assembles rigid inertia, hydrodynamic added mass and radiation, hydrostatic restoring, excitation, and linear PTO forces. The heave, fixed-hinge, and floating-joint layouts run the paired cases above. A fourth `linear_subspace` layout maps independent coordinates into arbitrary bodies; its RM3 two-body heave, Sphere free decay, and configured Sphere PTO checks agree with MATLAB. Arbitrary Simscape layouts, moorings, nonlinear hydro, and other application cases remain unsupported. |

The targeted [RM3 sea-state matrix export](https://github.com/cmudrc/wec-sim-python/actions/runs/37682612044)
confirms that MATLAB's applied added-mass matrix and adjusted rigid mass sum
to the original rigid mass plus HDF5 infinite-frequency added mass in the
active surge/heave/pitch block. The paired matrix gate passes for both bodies.
The split preserves the effective active mass exactly, so changing static
added-mass coefficients would have targeted the wrong mechanism.
The ordinary RM3 convolution solver retains implicit added mass. Paired
MATLAB tests also exercise `added_mass_scheme="simulink_delay"` for source
numerics; this compatibility setting is not the default dynamics path.
For Cases 3 and 4, paired gates now run the implicit default independently
of the explicit source-numerics option. The default has maximum float-surge
differences of 56.9 and 48.3 mm against MATLAB; the source option reduces
them to 28.0 and 24.8 mm. This comparison does not change the default
hydrodynamic coefficients or added-mass treatment.
The applied added-mass force in all three saved sea states uses acceleration
extrapolated from the two preceding 0.1 s samples to the current time minus
the Simulink block's `1e-7` s Transport Delay. After accounting for the
postprocessed pitch-inertia correction, that reconstruction matches the
reported force within `7e-7` N or N m across both bodies. This establishes
the source's acceleration feedback at output times. The optional
`simulink_delay` scheme applies the same mass split and delayed feedback.
The exact MATLAB joint adds `PTO stroke * sin(pitch)` to the two bodies'
relative surge, a term absent from the reduced Python geometry. In the saved
MATLAB traces this term peaks at 0.304, 0.222, and 0.140 m across the three
sea states, compared with 0.039 m in regular MCR case 8. Exact geometry alone
regressed regular-wave gates. Combining it with the source's delayed mass
feedback passes both the imported sea states and the previously paired
regular-wave cases. The delay is a Simulink numerical setting, not a measured
WEC property, and is confined to explicitly selected RM3 convolution
comparisons here.

The published RM3 body-to-body Cases 5 and 6 use the same suspect fitted
state-space radiation as the fourth `Radiation_Force_Options` setting. At 400 s,
MATLAB float surge is +0.515 m in Case 5 versus +0.110 m in convolution Case 3,
and +0.686 m in Case 6 versus +0.075 m in convolution Case 4. The fit's
effective common-surge damping is negative at low frequency although source
BEM damping is nonnegative at its sampled frequencies. Matching the MATLAB
state-space trajectory is a diagnostic reproduction, not evidence of physical
validity; it is excluded from validated parity until the fit is resolved.
The published Cases 5 and 6 change only the state-space radiation flag from
their respective convolution cases, apart from formatting of the inputs.
Over an 8 s window ending at 392 s, MATLAB float-surge means are 0.465 m
and 0.597 m for Cases 5 and 6, compared with 0.101 m and 0.047 m for
Cases 3 and 4. Repeating Cases 5 and 6 with a 0.05 s step gives means of
0.451 m and 0.592 m, so halving the published step does not remove the drift.
At zero frequency the fitted common-surge damping applied by MATLAB is
−14.1 kN s/m without cross-body coupling and −17.5 kN s/m with it.
The pinned MATLAB `bodyClass` sets the state-space direct term to zero even
though the BEMIO HDF5 stores nonzero terms. Including those saved terms in a
diagnostic transfer calculation reduces the negative zero-frequency values
to −3.01 and −3.45 kN s/m, but does not make the fits passive. This is a
low-frequency fit problem consistent with the drift, not proof that the
full nonlinear trajectory has only one cause.
Diagnostic PR #6 was closed without merging its state-space solver. The
convolution results for Cases 3 and 4 were identical, sample for sample,
between the original branch point and that diagnostic branch.
An independent no-wave check also separates the models: from a 0.01 m/s
common-surge perturbation, the default convolution solver ends at 0.00959 m/s
without cross-body radiation and 0.00934 m/s with it after 400 s. The pinned
fit's zero-wave linearization instead grows to 0.077 and 0.115 m/s, with
positive real growth rates of 0.0053 and 0.0064 s⁻¹. This is a regression
check on the default dynamics, not a Case 5–6 parity claim. The explicit
`simulink_delay` option is used only for source-numerics comparisons and does
not enable state-space radiation.

The paired case-runner checks compare the MATLAB and Python time grids, all
active body positions and velocities, and stationary degrees of freedom. RM3
checks also compare PTO stroke,
speed, force, and mechanical power; OSWEC checks compare PTO angle, angular
speed, torque, and mechanical power. MATLAB reports absorbed PTO power with a
negative sign, while the Python API reports positive absorbed power. The
Sphere checks now require a stationary 0 m case and narrower motion and force
limits for the displaced cases. Tolerances remain explicit in
`tests/test_case_dynamics_parity.py` and are calibrated against the pinned
MATLAB reference output, not against historical Python fixtures.

Source comparisons: [MATLAB wave class](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/waveClass.m),
[MATLAB wave-number function](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/functions/BEMIO/calcWaveNumber.m),
and [MATLAB body class](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/bodyClass.m).

## Reference model and case inventory

`tests/reference_cases.csv` inventories 63 published input files at the
pinned revisions: two core examples and 61 WEC-Sim Applications cases. It
identifies 21 RM3, 13 OSWEC, and 12 Sphere cases by explicit hydrodynamic
file paths, plus 17 other or dynamic cases. Generate the CSV with
`python tools/build_reference_case_inventory.py CORE_CHECKOUT APPLICATIONS_CHECKOUT`.
This is a source inventory, not a claim that every case runs in Python.

The `MATLAB reference model baselines` workflow runs the two canonical core
examples, all five Sphere free-decay cases, RM3 body-to-body Cases 1–4,
both RM3 PTO extension free decays, four RM3 radiation-force baselines (three
paired to Python),
all eight physical conditions from RM3 Multiple Condition Runs Option 1,
three actual `wecSimMCR` runs, each with eight saved body and PTO traces,
the published Sphere passive, reactive PI, declutching, and latching
controllers, a configured Sphere PTO case derived from the passive input,
and the published OSWEC `Nonhydro_Body` case. An earlier
[nine-job run](https://github.com/cmudrc/wec-sim-python/actions/runs/37641710801)
passed against fresh MATLAB R2025b outputs. The Option 1 MCR input was also
expanded into eight scalar simulations. The three driver baselines run
MATLAB's MCR orchestration and postprocessing directly. The pinned Option 2
workbook and Option 3 MAT file select the same ordered cases as Option 1.
Python's MCR interface pairs its case trajectories, PTO signals, average
powers, and power-matrix placement to each driver's actual output.
It records time, position, velocity, total force, and excitation force for
each body as CSV artifacts. New runs also record position, velocity, internal
mechanics force, and power for each PTO. The Applications source is pinned to
[`d53d4d4c9eda2581f04204f5d394a6ef84bb099e`](https://github.com/WEC-Sim/WEC-Sim_Applications/tree/d53d4d4c9eda2581f04204f5d394a6ef84bb099e).
The first run [passed all three model jobs](https://github.com/cmudrc/wec-sim-python/actions/runs/37471757955),
producing nine finite body trajectories with 4,001 samples each (two RM3,
two OSWEC, and five Sphere). The `1m-ME` and `1m` Sphere trajectories were
identical in the saved motion and force signals: the Morison element in that
published case has nonzero coefficients only in x while the free-decay motion
is in heave.
The renewed [model baseline run](https://github.com/cmudrc/wec-sim-python/actions/runs/37472556464)
also generated the Sphere HDF5 file with current MATLAB BEMIO and verified
Python preprocessing of it. `wecsim/linearHeave.py` uses that
preprocessing, the heave restoring coefficient, infinite-frequency added
mass, and radiation impulse-response kernel. It integrates the resulting
linear convolution equation with a fixed 0.01 s trapezoidal step. The five
Sphere trajectory comparisons [passed in the same MATLAB job](https://github.com/cmudrc/wec-sim-python/actions/runs/37476318082)
so that the source, HDF5, time grid, and output remain paired.
The [current reference-model run](https://github.com/cmudrc/wec-sim-python/actions/runs/37483219008)
passed all RM3, OSWEC, and Sphere jobs, including Python preprocessing of
the current RM3 and OSWEC HDF5 inputs for both bodies and the focused dynamics
comparisons described above.
The RM3 baseline also checks a reduced Python heave model and a coupled
surge/heave/pitch model. The OSWEC baseline sets the published case's
`waves.phaseSeed` to 1 so repeat runs use one reproducible realization,
saves its realized phase matrix, and checks Python PM binning, wave elevation,
directional excitation, and hinged-pitch motion against that same realization.
The Python wave generator can use its own integer seed for standalone runs,
but its random sequence differs from MATLAB's Threefry generator. The
case-driven runner shares a generalized dynamics engine across the validated
layouts. A PM phase CSV can replay the actual MATLAB realization.
The RM3 body-to-body job generates the Applications HDF5 input with the pinned
MATLAB BEMIO code, then compares the Python solver with coupling off and on
against paired MATLAB Cases 1 and 2. The CLI selects these hydrodynamic modes
with `rm3` and `rm3 --b2b`; neither command executes an arbitrary application
input file.
The [four-job reference-model run](https://github.com/cmudrc/wec-sim-python/actions/runs/37486493413)
passed RM3, RM3 body-to-body, OSWEC, and Sphere, including the Sphere CLI
smoke test against the MATLAB-generated HDF5 input.
The [case-driven dynamics run](https://github.com/cmudrc/wec-sim-python/actions/runs/37491350883)
passed all four jobs again using JSON cases through the main runner. It also
compared the mapped linear-coordinate RM3 heave and Sphere free-decay cases
against the paired MATLAB trajectories.
Optional linear PTO equilibrium offsets and scalar pretension use the same
generalized dynamics engine. Their nonzero-force behavior has analytical and
case-level tests; the paired MATLAB reference cases use zero offsets.
The configurable small-motion device path also accepts named body motions and
body-local or fixed-world PTO endpoints. Its attachment geometry, projected
stroke, generalized forces, and damping power have analytical and case-level
checks. The configured Sphere case pairs a nonzero spring and extra damper
with MATLAB output while specifying a shifted body-local PTO point. Because
that body only heaves, the shifted point does not affect its trajectory;
attachment-location dynamics remain unpaired.
`wecsim.WEC` provides a Python builder for this same validated path and
returns named NumPy body, coordinate, and PTO histories. The JSON case runner
remains available for saved cases; the Python builder currently covers the
`linear_subspace` layout only.
`python -m wecsim CASE.json --output motion.csv` runs a
supported dynamics configuration. The case declares wave, body, constraint,
PTO, and time settings; the result includes all six body position and
velocity coordinates, applicable wave and PTO signals, and a JSON record of
the case and HDF5 hashes, NumPy version, and Git state. The older
`referenceRunner` remains as a preset CLI for the three model families.
Neither command parses an arbitrary WEC-Sim or Simscape input file.
Expanded application cases remain inventoried. The manually dispatched
`MATLAB reference application regression`
workflow can run the upstream test suites for all 18 application folders
containing the 46 explicit RM3, OSWEC, and Sphere cases. This checks the
pinned MATLAB reference for those cases, but there is no Python time-series
comparison yet.

The [first full application sweep](https://github.com/cmudrc/wec-sim-python/actions/runs/37472556454)
recorded 42 passed methods, two failed Passive Yaw assertions, two MoorDyn
methods filtered by upstream's CI assumption, and an empty Multiple Wave
Spectra suite. Both Variable Hydro methods passed under MATLAB R2025b. A
targeted direct run subsequently passed Multiple Wave Spectra, bringing
the exercised passing checks to 43.

The first application sweep revealed upstream test-suite limitations. In
`Passive_Yaw`, two stored irregular-yaw regression checks fail; the same two
checks fail in [the Applications repository's own CI run](https://github.com/WEC-Sim/WEC-Sim_Applications/actions/runs/36746373598)
on both Linux and Windows under MATLAB R2024b. MoorDyn tests in `Mooring` and
`Paraview_Visualization` explicitly skip themselves on GitHub CI. Forcing the
MoorDyn test exposed a native library mismatch (`GLIBCXX_3.4.32` unavailable
from MATLAB's bundled `libstdc++`), so those cases remain unverified in CI.
The `Multiple_Wave_Spectra` test class is excluded by MATLAB because its
class name does not match its filename. Our harness generates its OSWEC HDF5
file with BEMIO, runs the input file directly, and requires body output
instead of accepting an empty test suite as success. The
[direct case passed](https://github.com/cmudrc/wec-sim-python/actions/runs/37476616198).
The Passive Yaw and MoorDyn findings are MATLAB-reference gaps, not evidence
of Python dynamics agreement or disagreement.

## Known differences and next reference case

Current MATLAB WEC-Sim has changed its PM/JS spectra, seeded phase generator,
object properties, and some wave inputs since the Python port was written.
The historical fixtures therefore cannot establish parity for those modes.
The focused OSWEC PM implementation above is checked against current MATLAB;
the inherited general `WaveClass` still needs corresponding updates for other
irregular-wave cases.
RM3 state-space Cases 5 and 6 remain a radiation-fit investigation: a
passive fit would need its own evidence against source BEM data and
convolution trajectories before it could count as physically validated.
The [pinned Applications `End_Stops` input](https://github.com/WEC-Sim/WEC-Sim_Applications/blob/d53d4d4c9eda2581f04204f5d394a6ef84bb099e/End_Stops/wecSimInputFile.m) sets
`upperLimitTransitionRegion` and `lowerLimitTransitionRegion` to 0.5 m, but
the [pinned `ptoClass.hardStops` defaults](https://github.com/WEC-Sim/WEC-Sim/blob/0753b2e47f2457c078751dcfe5d251d1767b80ab/source/objects/ptoClass.m) and the referenced translational PTO
Simulink block read the fields ending in `TransitionRegionWidth`. The input
therefore leaves the effective widths at their 1e-4 m defaults. The
[fresh source export](https://github.com/cmudrc/wec-sim-python/actions/runs/37689167739)
confirms these effective settings and supplies body/PTO trajectories,
hydrodynamic force components, and original/applied mass matrices. The PTO
stroke reaches −0.684 to +0.687 m versus about ±0.865 m in the published
no-stop case, so contact materially changes the trajectory. A refined 400 s
source run records two samples inside the effective transition width; the
Python `LinearHardStops` smoothing law reconstructs its logged PTO force to
within `2.4e-7` N. The upstream
[test](https://github.com/WEC-Sim/WEC-Sim_Applications/blob/d53d4d4c9eda2581f04204f5d394a6ef84bb099e/End_Stops/TestEndStops.m)
only checks that `wecSim` runs.

The source regular-wave case applies split added mass with delayed
acceleration feedback. Reconstructing its translational added-mass force from
the two preceding 0.1 s output accelerations agrees to numerical precision
until 84.4 s; 78 later output samples differ by more than 0.001 N, with a
maximum residual of 17.9 MN. Those intervals are consistent with internal
variable solver steps near impacts, whose acceleration history is absent from
the exported 0.1 s samples. The [0.05 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37691660845)
and [0.025 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37692499385)
and [0.0125 s source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37693255891)
have no such saved-step residual above 0.001 N through 120 s. PTO stroke
changes by up to 20.7 mm from 0.05 to 0.025 s and 2.79 mm from 0.025 to
0.0125 s, versus 149.4 mm from 0.1 to 0.05 s. Ordinary PTO damper energy
over 120 s falls from 11.176 MJ in the published 0.1 s run to 10.518 MJ at
0.0125 s.

Over the full 400 s, the two refined MATLAB runs differ by at most 2.79 mm
PTO stroke, 18.8 mm/s PTO speed, and 242.4 kN PTO force. Their ordinary
damper energies are 55.062 and 54.823 MJ, a 0.437% difference. Neither
refined run has a saved-step added-mass recurrence residual over 0.001 N;
the published 0.1 s run has 78 such samples and absorbs 60.941 MJ. The
source convergence gate also checks both bodies' active positions and
velocities. These are numerical step checks, not altered WEC settings.

The Python hard-stop model instead applies added mass implicitly and resolves
contact forces between output samples. Its 0.025 and 0.0125 s output settings
change PTO stroke by less than 1 µm over 400 s. This physical dynamics
path is available through the RM3 solver and case runner; it does not embed
the source's delayed acceleration loop. Against the 0.0125 s MATLAB run, the
Python body's largest surge and heave position errors are 1.08 and 0.936 mm;
the largest pitch error is `6.74e-5` rad. Ordinary damper energy is 54.747 MJ,
0.0762 MJ (0.139%) below the refined source. General nonlinear stop damping
and arbitrary WEC layouts remain outside this RM3 validation.
Further dynamics targets include other published OSWEC and Sphere
configurations. The Sphere free-decay cases already have direct Python motion
comparisons. A comparison must
record both code revisions, the HDF5 input, time step, outputs, and numerical
tolerances.

The inherited `paraviewClass.py` was invalid Python mixed with unfinished
MATLAB code. The wave-surface calculation is now importable; VTP serialization
and body visualization explicitly raise `NotImplementedError` until ported.
The older object tests import duplicate copies of classes inside test folders.
The parity tests here import the production files instead.
Six inherited files under `tests/test_simulink` still contain unfinished
MATLAB-like Python and do not parse; they are outside the collected pytest
suite. The `wecsim` package compiles, and the scoped production parity suite runs.
