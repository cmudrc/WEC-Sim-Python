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
| Current PM and JONSWAP spectra and replayed elevation | Pinned MATLAB `waveClass` with both Traditional and EqualEnergy 64-frequency discretizations, 8 s peak period, significant heights of 2.5 m and 4 m, and its Threefry phase output; [live R2025b gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37777444325) | Production `WaveClass.waveSetup` uses the current height-dependent PM and JONSWAP equations, including inferred JONSWAP `gamma`. All six live checks pass. In the new equal-energy pairs, the maximum difference across frequency, width, and spectrum entries is `5.4e-15`; deep-water power differs by at most `2.2e-11` W/m, and origin and three-marker elevation by at most `7.4e-15` m. PM uses three incident directions. Python replays the exported phase matrix because NumPy and MATLAB do not generate identical phases from the same integer seed. |
| Narrowed PM and JONSWAP frequency ranges | Same pinned MATLAB `waveClass`, with `bem.range=[0.5,1.5]` rad/s inside a `[0.4,2.0]` rad/s BEM interval; [passing R2025b gate](https://github.com/cmudrc/wec-sim-python/actions/runs/37780011107) | The public `PMWave` and `JONSWAPWave` settings now pass an optional range through the hydrodynamic `WEC` runner. The production component builders' selected frequencies, widths, and spectral values differ from MATLAB by at most `5.8e-15`; replayed origin elevation differs by at most `1.9e-14` m. The live gate also checks power and three marker elevations through `WaveClass`. Out-of-BEM range endpoints clamp to HDF5 limits as in MATLAB; hydrodynamic bodies retain HDF5 water depth. This validates wave preprocessing and API propagation, not a new coupled-motion case. |
| Irregular Bretschneider equal-energy binning | Original MATLAB-era test constants | Production `WaveClass` agrees after NumPy/SciPy API updates. This is historical compatibility only: current MATLAB WEC-Sim rejects the `BS` option. |
| Wave-surface grid for no-wave, regular, and directional irregular waves | Current MATLAB `waveClass.waveElevationGrid` equations | Implemented and checked for regular and two-direction irregular cases. |
| RM3 HDF5 hydrodynamic input | Original `rm3.h5` | Both bodies load, including their names and water depth, under NumPy 2. |
| Current RM3 and OSWEC HDF5 hydrodynamic inputs | Pinned current MATLAB core examples | Both bodies in each model load and complete regularCIC or directional irregular preprocessing with finite restoring, added-mass, excitation, and radiation arrays. This is input compatibility, not motion parity. |
| RM3 regular-wave force preprocessing | Original MATLAB-generated `body_1_test` constants and files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF; the RM3 runner now completes preprocessing. |
| RM3 irregular-wave force preprocessing | Original MATLAB-generated `body_2_test` files | Production `BodyClass` agrees for restoring stiffness, added mass, excitation, and radiation IRF. |
| RM3 body interaction force preprocessing | Original MATLAB-generated `body_4_test` through `body_9_test` files | All six regular/regularCIC variants, including body-to-body coupling on/off and state-space radiation on/off, agree for both RM3 bodies on restoring stiffness, added mass, excitation, radiation IRF, and state-space matrices where present. |
| OSWEC directional irregular force preprocessing | Original MATLAB-generated `body_3_test` files | Production `BodyClass` agrees for restoring stiffness, added mass, radiation IRF, and three-direction excitation after replacing removed SciPy `interp2d`. |
| Sphere `noWaveCIC` heave free decay (0 m, 1 m, 1 m with Morison elements, 3 m, 5 m) | Current MATLAB WEC-Sim and MATLAB-generated Sphere HDF5 | The focused Python linear heave solver agrees over 40 s to maximum differences of 0.45 mm position, 0.64 mm/s velocity, and 363 N total force in the 5 m case. The Morison element has only x-direction coefficients in the published 1 m case, so it does not affect heave. |
| Sphere sampled-elevation heave | Pinned `Free_Decay/0m` Sphere model with only its no-wave input replaced by a reproducible two-frequency `elevationImport` record; [fresh R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37785332545) | The public `WEC.run(ImportedElevationWave(...))` path reuses the production imported-elevation convolution and runs all 4,001 samples over 40 s. With the source body block's second force ramp explicitly selected, maximum MATLAB/Python differences are `6.7e-15` m elevation, `3.7e-9` N or N m across six excitation components, 0.0243 mm heave, and 0.0323 mm/s speed. Paired gates are `1e-12` m, `1e-6` N or N m, 0.1 mm, and 0.1 mm/s. This derived case validates one-body heave dynamics with a sampled sea; the published two-body RM3 MooringMatrix imported-elevation motion remains separately paired through its floating-joint solver. |
| Fixed monopile Cartesian Morison force | Pinned MATLAB Applications `Morison_Element/morisonElement` and R2025b run, with no HDF5 bodies | The public `WEC.fixed_body` and `WEC.morison_element` configuration replays the published 400 s, three-heading PM sea and its six-component stationary-body force. All 500 equal-energy frequencies, spectral amplitudes, widths, and finite-depth wavenumbers agree within `1e-12`; elevation agrees within `1.9e-13` m. On all 40,001 samples, the largest force/moment difference is `7.3e-6` N m against a `1e-3` paired gate, and both body positions and velocities agree exactly. The source logs the negative of physical Morison force. Its `irregWaveMorison.m` uses the first random-phase column for every force heading even though wave elevation uses heading-specific phases; Python's explicit `phase_mode="matlab_shared"` reproduces that source behavior while the default uses heading-specific phases. This validates the published fixed Cartesian element, not moving elements or normal/tangential coefficient mode. |
| Fixed hydrodynamic monopile and tower | Pinned MATLAB Applications `Morison_Element/monopile`, generated `monopile.h5`, and [fresh R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37858193863) | The public `run_fixed_hydro_monopile` synthesizes the published 500-bin PM sea, six-component hydrodynamic excitation, two stationary-body force traces, and two fixed-joint reactions. On all 40,001 samples over 400 s, maximum elevation error is `2.1e-13` m, excitation and hydro-body total-force error `5.9e-7` N or N m, and joint-reaction error `2.6e-7` N or N m, against paired gates of `1e-10` m and `1e-5` N or N m. Both positions and velocities agree exactly. Source radiation, added mass, Morison force, damping, and acceleration are zero in this stationary case. The tower's weight and joint reactions close static force and moment balance. This verifies the published fixed hydrodynamic case, not a moving monopile or general constraint solver. |
| MBARI `Cable` coupled three-body regular-wave case | Pinned Applications `Cable`, generated `mbari.h5`, and [fresh R2025b force and joint export](https://github.com/cmudrc/wec-sim-python/actions/runs/37859525482) | The public `run_mbari_cable` independently advances both hydrodynamic bodies and the nonhydrodynamic cylinder with the published spherical joint, axial spring/damper, and drag on both cable endpoint bodies. Over all 12,001 samples and 120 s, maximum body-position error is 4.3 cm, pitch error 0.0084 rad, and body-speed error 0.422 m/s at a stiff cable snap. Cable force differs by 3.03 kN RMS and 39.9 kN at its worst instant against a 200 kN source peak; signed force-impulse error is 234 N s. The paired gates require positions below 5 cm, active speeds below channel-specific limits up to 0.45 m/s, cable force below 45 kN pointwise and 3.5 kN RMS, and impulse error below 500 N s. The source radiation, hydrostatic, and viscous body force laws reconstruct on saved states within `1e-7` N or N m; the source HDF5 added-mass and damping matrices agree with Python within `1e-8`. Given the saved body poses, rates, and accelerations, cable axial force, endpoint drag, and each 1 kg endpoint inertia reconstruct the source lower/upper vertical joint loads within 2/7 N. The independent runner retains physical implicit added mass and omits the two 1 kg endpoint inertias; the pinned Simscape source shifts added mass into rigid bodies with delayed acceleration feedback. This is planar published-case trajectory parity, not a general six-DOF cable or constraint solver. |
| WaveBot `Load_Mitigating_Controls/CalcImpedance` three-DOF system identification | Pinned Applications multisine A, WAMIT/BEMIO `waveBotBuoy.h5`, and a fresh R2025b run of the published 600 s no-wave case | `run_wavebot_impedance` independently integrates surge, heave, and pitch with 10 s radiation memory, implicit added mass, mooring stiffness/damping, quadratic drag, and the published multisine force. Source-convention motion differs from MATLAB by at most 0.13 mm surge, 0.20 mm heave, and 0.00028 rad pitch over all 60,001 samples; speed differences stay below 1.7 mm/s. Source force diagnostics reconstruct actuation, mooring, quadratic drag, and linear damping; the independent mooring-force error stays below 3.1 N. The published input uses `linearDamping(1:2:5)`, which fills the first matrix column in MATLAB and couples heave/pitch damping to surge speed. Its pitch command is applied with a negative sign. `source_linear_damping=True` reproduces that source input; the default uses a dissipative diagonal matrix and differs from the MATLAB heave/pitch trajectories by up to 17.7 mm/0.0079 rad. The source damping matrix delivers positive instantaneous power during 42% of samples, although net work is dissipative over the run. This covers open-loop system identification, not the separate adaptive `ControlTests` case. |
| Fixed monopile with uniform, power-law, and linear current | Three derived 20 s, one-heading cases from pinned `Morison_Element/morisonElement`; [R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37820035413) | `PMWave(current=Current(.8, 45, profile, 30))` adds a ramped horizontal current to the fixed element's fluid velocity. Each profile has a full 2,001-sample MATLAB wave and six-component force time-series pair. The largest force/moment errors are `1.23e-6`, `1.23e-6`, and `1.16e-6` N or N m for uniform, one-seventh power, and linear profiles, respectively, against a `1e-3` gate. Wave elevation passes a `1e-10` m gate for each. Active current with multiple wave headings is rejected because the pinned MATLAB force function counts it once per heading. Moving-body currents remain unsupported. |
| Moving Cartesian Morison source force in regular waves | Pinned `regWaveMorison.m` option 1, called at eight prescribed six-DOF states during the [MORISON_FIXED MATLAB baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37768083149) | The separate `regular_morison_source_force` diagnostic matches all six MATLAB force and moment components within `5.5e-12` N or N m, including nonzero body velocity and acceleration, angular motion, wave ramp, and an emerged zero-force state. These are prescribed-state force checks, **not coupled moving-body trajectory parity**. The pinned source rotation matrix is nonorthogonal at the tested nonzero attitudes and its angular kinematics use the unrotated local point; the diagnostic reproduces those source rules without changing the default WEC dynamics. General moving Morison layouts remain unsupported in the public `WEC` runner. |
| Sphere moving Morison heave free decay | The pinned `Free_Decay/1m-ME` application with its originally surge-only element changed to axial heave coefficients (`Cd=Ca=1`, area 100 m², volume 20 m³); [fresh R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37770321922) | The public `WEC.morison_element` couples body-local axial drag and acceleration-dependent added mass to a hydrodynamic heave body in still water. Over all 4,001 samples, independent Python motion differs by at most 0.325 mm heave and 0.952 mm/s speed. After the first second, physical Morison force differs by at most 43.2 N; the 40 s signed force-impulse difference is 198 N s. MATLAB logs its acceleration-feedback Morison force as zero at the first two samples and has a 41.2 kN pointwise difference from Python during startup; reconstructing MATLAB force from its *saved* acceleration only agrees within 14 N after 0.5 s. Python keeps implicit added mass rather than reproducing that source feedback transient. This derived case validates coupled single-body heave in still water. |
| Sphere moving Morison heave in regular waves | The same pinned `Free_Decay/1m-ME` Sphere case, with heave coefficients above and a 1 m, 8 s regular wave ramped over 10 s; [passing R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37773024911) | The public single-body heave solver evaluates fluid velocity and acceleration at the moving body-local element, uses relative velocity for drag, and includes the element's added mass in the acceleration solve. Across 4,001 samples over 40 s, wave elevation differs by at most `3.8e-15` m, heave by 2.26 mm, and speed by 2.52 mm/s. After the first second, physical Morison force differs from the logged source force by at most 0.935 kN over a roughly 65 kN peak. The source floating joint also permits surge and pitch: its surge reaches 1.42 m and pitch 0.0293 rad, whereas the Python configuration only moves in heave. On the source's actual six-DOF states, direct `regWaveMorison.m` replay matches the logged force within 10.1 N after startup; projection onto heave alone changes that source force by up to 0.973 kN. The Python heave force matches the pinned source function on those projected states within `5.9e-10` N, against a `1e-6` N gate. This is **reduced heave comparison**, not full three-DOF trajectory parity. The initial source acceleration-feedback force spike is excluded from the post-startup force gate. Other moving layouts and currents remain unsupported. |
| Sphere moving Morison in three-DOF regular waves | The derived 1 m/8 s `Free_Decay/1m-ME` case with its published surge/heave/pitch floating joint; [fresh passing R2025b run](https://github.com/cmudrc/wec-sim-python/actions/runs/37776059648) | The public `WEC` integrates a body-local axial element through all three active coordinates using a proper pitch rotation, relative fluid velocity, fluid acceleration, and an implicit positive-semidefinite added-mass matrix. Against the 4,001-sample, 40 s MATLAB trajectory, the independent Python run differs by at most 12.2 mm surge, 1.01 mm heave, and `9.8e-5` rad pitch; velocities differ by at most 1.32 mm/s surge and heave, and `4.6e-5` rad/s pitch. After the first second, force and moment maxima differ by 429 N surge, 149 N heave, and 176 N m pitch, with explicit paired gates on each channel. Evaluating the physical Python force law on the saved MATLAB states separates source rotation conventions from trajectory error: differences are at most 427 N surge, 91 N heave, and 180 N m pitch. The source-specific `regWaveMorison.m` diagnostic remains separate from the physical default. This validates the published three-coordinate Sphere case; general six-DOF moving elements, currents, and other coefficient layouts remain unsupported. |
| Sphere `Controls/MPC` closed-loop heave and PTO | Pinned MATLAB Applications `Controls/MPC` at `d53d4d4`, pinned core at `0753b2e`, and [instrumented R2025b baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37743807105) | The public `run_sphere_mpc` solves the 400 s published case independently with 500 JONSWAP bins, a fourth-order radiation fit, 201-row AR excitation forecast, constrained 31-variable QP, delayed force-rate actuation, and physical heave radiation memory. Replayed phase frequencies, spectra, and widths agree within `1e-12`, wave elevation within `3e-12` m, and heave excitation within `1e-7` N. All prediction and QP matrices agree within `1e-10`; all 597 forecasts after a full history differ by at most `5e-7` N, and all 391 active QP first commands agree within 1 N/s using MATLAB's logged state and forecast. In an independent Python closed-loop run, maximum 40,001-sample differences are 1.011 mm physical heave, 1.501 mm/s speed, 17.2 N/s command rate, 16.8 N PTO force, and 2.89 kW source-signed instantaneous power. The MATLAB controller's internal plant is distinct from the physical body; its position is relative to equilibrium. The 2.5 m/8 s input infers JONSWAP `gamma=1`. The source's actual PTO force reaches 2.39 MN despite its 2 MN predicted constraint. This gate covers this single-body, heave-only application and published settings; general multi-body MPC remains unsupported. |
| Sphere `Mean_Drift` regularCIC application | Pinned MATLAB Applications `Mean_Drift` input and its BEMIO-generated Sphere HDF5 | The public Python `WEC` builder selects the control-surface mean-drift coefficient and runs the published 100 s surge, heave, and pitch case with 10 s radiation memory. The source excitation equals the first-order regular force plus amplitude-squared mean drift, each with one ramp factor; the surge drift component is 12.247 N after ramp. Paired maximum position differences are 0.043 mm surge, 0.551 mm heave, and `3.2e-7` rad pitch; velocity differences are below 0.081 mm/s surge, 1.738 mm/s heave, and `1.7e-7` rad/s pitch. The wave elevation matches within `3e-15` m, and reconstructed radiation forces differ by at most 0.24 N surge, 2.95 N heave, and `0.0002` N m pitch. The MATLAB linear model moves 16.41 m in surge by 100 s, so agreement is numerical parity, not evidence that the small-motion hydrodynamics remain physically accurate over that travel. |
| RM3 regular-wave heave, excitation, and PTO | Current MATLAB RM3 example and current RM3 HDF5 | A focused two-body linear heave solver with relative-motion PTO damping agrees locally over 400 s within 5.9 mm position and 6.1 mm/s velocity for both bodies. Heave excitation matches to less than `1e-5` N; PTO internal force differs by at most 8.2 kN over a 1.63 MN range. This reduced model does not cover surge, pitch, or full Simscape joint mechanics. |
| RM3 coupled surge, heave, pitch, and PTO | Current MATLAB RM3 example and current RM3 HDF5 | A four-coordinate two-body model uses the published pitched-slider geometry, body inertias, regular-wave forcing, and nonlinear rotation kinematics. Over 400 s, the two body surge positions differ by at most 36.8 mm, heaves by 4.18 mm, shared pitch by 0.000672 rad, and PTO force by 4.74 kN. It covers the canonical active DOFs but not general Simscape joint mechanics or other RM3 cases. |
| RM3 body-to-body Cases 1 and 2 | Pinned MATLAB Applications inputs and RM3 HDF5 generated with current BEMIO | Both published regular-wave cases run for 400 s with coupling off and on. Across the two, body positions differ by at most 36.8 mm surge, 4.57 mm heave, and 0.000948 rad pitch; PTO force differs by at most 5.38 kN over a 1.65 MN range. Turning on cross-body coupling substantially reduces the body-1 heave error against Case 2. These are reduced-model comparisons, not full Simscape mechanics. |
| RM3 body-to-body Cases 3 and 4 | Pinned MATLAB Applications `regularCIC` cases, with coupling off and on | Both the ordinary implicit-mass solver and the opt-in `simulink_delay` comparison are paired over 400 s. The ordinary solver differs by at most 67.2 mm surge, 3.64 mm heave, 0.00162 rad pitch, and 4.78 kN PTO force across both bodies and cases. The source-numerics setting reduces those maxima to 28.2 mm, 0.253 mm, 0.0000325 rad, and 0.292 kN; PTO stroke, speed, and power then differ by at most 0.235 mm, 0.244 mm/s, and 0.392 kW. Cases 5 and 6 use MATLAB's fitted radiation state-space model and are not validated by this convolution solver. |
| RM3 MooringMatrix with imported elevation | Pinned MATLAB Applications `Mooring/MooringMatrix` input, its 3,600 s `etaData.mat`, and BEMIO-generated RM3 HDF5 | The public `run_case` interface loads the MAT elevation, applies the configured joint surge spring, and runs the full 400 s comparison. The source wave agrees with Python interpolation and ramp to `8e-14` m. Production `BodyClass` uses the source's linear excitation-IRF interpolation and noncausal `same` convolution; an explicit opt-in reproduces the source body block's additional force ramp. All logged excitation-force components agree within `1e-3` N or N m. The logged mooring force equals `−100,000` N/m times its surge displacement within `1e-5` N. On all 40,001 samples, implicit-added-mass dynamics agree within 9.77 mm surge, 2.53 mm heave, 0.000605 rad pitch, 2.55 mm PTO stroke, 0.279 mm/s PTO speed, 334 N PTO force, and 1.17 kN mooring force. General mooring layouts remain unsupported. |
| RM3 End_Stops force law and refined-step trajectory | Pinned MATLAB Applications `End_Stops` input; R2025b runs at published 0.1 s and diagnostic 0.025/0.0125 s steps through the published 400 s duration | Effective stroke bounds are ±0.6 m, each spring is 100 MN/m, stop damping is zero, and the 0.0001 m transition is active. The source's logged stop force matches the Python smooth-stop law within `1e-5` N even at saved samples inside the transition. The 0.025-to-0.0125 s source stroke difference is at most 2.79 mm over 400 s; the published 0.1 s run differs from the refined source by up to 286 mm after contact. The opt-in Python adaptive solver with implicit added mass agrees with the 0.0125 s MATLAB trajectory over 400 s within 0.986 mm PTO stroke, 6.77 mm/s PTO speed, 98.2 kN PTO force, and 0.139% ordinary damper energy. Paired gates cover both active body positions and velocities, PTO motion and force, and integrated energy; a separate source convergence gate checks both refined MATLAB steps. |
| RM3 radiation force options: constant, convolution, and FIR | Pinned MATLAB Applications `Radiation_Force_Options`, BEMIO-generated RM3 HDF5, 12 s regular waves | All three published 500 s settings are checked against the Python floating-joint solver. The FIR path uses sampled kernel taps with a held radiation force during each RK4 step. Across the three cases, the largest differences are 29.6 mm surge, 2.11 mm heave, 0.00140 rad pitch, 3.25 mm PTO stroke, 1.55 mm/s PTO speed, 1.86 kN PTO force, and 1.57 kW PTO power. Direct FIR radiation-force differences remain under their paired gates. MATLAB FIR differs from MATLAB convolution by up to 345 mm surge, so the methods are not interchangeable for this case. The application's state-space setting is preserved only in closed, unmerged diagnostic PR #6 because of negative low-frequency damping in its fit. |
| RM3 PTO extension float and spar free decays | Pinned MATLAB Applications `RM3_PTO_Extension` inputs and BEMIO-generated RM3 HDF5 | The no-wave floating-joint solver initializes the float at +5 m or the spar at −5 m, reproducing each published case's +5 m PTO stroke. Over 30 s, maximum active-body heave differences are 47.2 mm for float and 6.85 mm for spar; velocity differences are 86.1 mm/s and 1.87 mm/s. Both initial poses and the passive body's heave agree within numerical precision; PTO force and power are zero in Python and within `1.7e-8` in MATLAB output. |
| RM3 Multiple Condition Runs Option 1 physical sweep | Pinned MATLAB Applications RM3 input, expanded to eight scalar height × period × PTO-damping combinations | The Python `regularCIC` floating-joint solver agrees across all eight 400 s conditions within 75.1 mm surge, 0.387 mm heave, 0.000183 rad pitch, 0.413 mm PTO stroke, 0.231 mm/s PTO speed, 0.479 kN PTO force, and 0.487 kW PTO power. These scalar runs validate each physical condition; all three MCR drivers are separately checked below. |
| RM3 Multiple Condition Runs Options 1–3 orchestration | Pinned MATLAB Applications arrays, Option 2 Excel grid, Option 3 MAT file, and three actual MATLAB R2025b `wecSimMCR` runs | The Python loaders produce the same ordered eight-case table for all three published inputs. Each driver has a separate paired gate for every body's 4,001 position and velocity samples, every PTO stroke, speed, force, and power trace, the eight mean powers, and two 2×2 power matrices with period rows and height columns. Python mean absorbed powers differ from MATLAB by at most 0.155 kW or 0.059% across eight conditions in the MAT driver. The Option 2 Excel and Option 3 MAT driver outputs agree numerically exactly; Option 1 array and Option 3 MAT body positions differ by at most `1.9e-12` m. MATLAB's PTO power sign is opposite Python's positive absorbed-power convention. Phase-seed sweeps, multiple PTOs, and other postprocessing are not yet paired. |
| Public RM3 floating-joint configuration | Pinned MATLAB MCR Option 1 array condition 1 and published MooringMatrix application; fresh paired [MCR](https://github.com/cmudrc/wec-sim-python/actions/runs/37787926107) and [mooring](https://github.com/cmudrc/wec-sim-python/actions/runs/37787936962) runs | `WEC.floating_joint` routes a Python-configured two-body device through the already paired pitched-slider solver and returns named body, four-coordinate, and PTO histories. In the 400 s MCR case, maximum differences are 35.1 mm surge, 0.228 mm heave, `3.46e-6` rad pitch, 0.243 mm PTO stroke, 163 N PTO force, and 118 W absorbed power; explicit gates cover positions, speeds, stroke, force, and power at all 4,001 samples. A second public configuration reproduces the 400 s imported-elevation MooringMatrix run, including its initial spar offset and joint surge spring, under the existing 40,001-sample body/PTO/mooring gates. Both fresh MATLAB jobs pass. The builder's implicit added-mass default remains physical; source-compatible delayed feedback is explicit for the MCR comparison. This layout does not represent arbitrary Simscape joints or off-axis PTO endpoints. |
| RM3 imported-spectrum sea-state MCR | Pinned MATLAB Applications `RM3_MCROPT3_SeaState` and three actual MATLAB MCR runs | All three imported component tables agree to `4.9e-15`, wave elevations to `1.0e-13` m, and active excitation forces and moments to `3.8e-8` N or N m. The public `WEC.run(ImportedSpectrumWave(...))` path now replays all 4,001 incident-wave samples in each 400 s case within `1.0e-13` m; its general WEC dynamics are not claimed to reproduce the published four-coordinate RM3 motion. On MATLAB body velocities, the 60 s radiation convolution matches every logged force component within `5.4e-9` N or N m. The applied added-mass force reconstructs from the exported matrix and delayed acceleration within `7e-7` N or N m. The exact pitched-slider geometry and source-compatible mass feedback bring all three 400 s paired trajectories within 14.4 mm surge, 1.03 mm heave, and 0.000638 rad pitch; velocity differences are below 3.26 mm/s surge, 0.883 mm/s heave, and 0.000256 rad/s pitch. PTO stroke, speed, force, and absorbed power differ by at most 2.44 mm, 1.41 mm/s, 1.69 kN, and 1.19 kW. Mean absorbed power differs by at most 29 W. Explicit paired gates passed for the specialized MCR solver; PR #18 was merged. |
| RM3 PTO-Sim direct linear generator | Pinned MATLAB Applications `PTO-Sim/RM3/RM3_DD_PTO`, its BEMIO-generated RM3 HDF5, and a fresh 400 s `ode4` run at 0.0005 s | Both the focused Python two-heave runner and a public `WEC.pto(linear_generator=...)` configuration integrate the source block's d/q flux states and electrical angle with float and spar motion. The public configuration uses the existing body-local endpoint and coordinate maps. Source-only checks verify PTO relative motion, friction, load voltage, and mechanical/electrical power identities. Across all 40,001 saved samples, both paired gates require body heaves/velocities and PTO stroke/speed within 10 µm or µm/s, generator force and powers within 0.005 N or W, and all three phase currents/voltages within 0.0001 A/0.01 V. The [fresh source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37749845161) and local paired tests pass. This validates the published vertical two-body regular-wave layout, not general PTO-Sim electrical networks, radiation-memory coupling, or arbitrary attachment geometry. |
| RM3 PTO-Sim rectified hydraulic PTO | Pinned MATLAB Applications `PTO-Sim/RM3/RM3_cHydraulic_PTO`, PTO-Sim library, BEMIO RM3 HDF5, and [fresh 400 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37823137925) | A dedicated Python two-heave runner couples the published float and spar, 60 s radiation IRF, PM sea, compressible cylinder, rectifying valve, gas accumulators, hydraulic motor, generator, and PI load controller. Replayed sea elevation agrees within `1e-10` m; source-driven component equations agree to numerical precision. Across 40,001 independently integrated samples, maximum float/spar heave errors are 2.59/4.31 mm, PTO stroke 4.87 mm, PTO speed 3.55 mm/s, cylinder force 14.64 kN against a 646 kN source peak, cylinder pressures 0.335 MPa, accumulator pressure 1.46 kPa, shaft speed 0.045 rpm, current 0.023 A, and voltage 0.023 V. Valve-port flow differs by up to 0.0442 m³/s at sharp reversals (mean absolute error 0.000856 m³/s); integrated port-volume error stays below 0.000434 m³. The paired gates cover these channels. This validates the published vertical-slider case, not general hydraulic layouts or attachment geometry. The source PI law commands negative load resistance after 95.35 s for 76.2% of samples; numerical agreement does not establish passive electrical operation. |
| Sphere passive controller | Pinned MATLAB Applications `Controls/Passive (P)` and MATLAB-generated Sphere HDF5 | The Python heave model represents the published proportional controller as an 860,870 N s/m damper. Maximum differences over 400 s are 0.063 mm position, 0.065 mm/s velocity, 56 N controller force, and 32 W controller power after accounting for MATLAB's opposite force and power signs. |
| Sphere reactive PI controller | Pinned MATLAB Applications `Controls/Reactive (PI)` and MATLAB-generated Sphere HDF5 | The published active feedback maps to signed linear PTO stiffness and damping in the existing Python heave coordinate. Against fresh MATLAB output over 200 s, maximum differences are 8.45 mm heave, 5.49 mm/s velocity, 4.85 kN controller force, and 35.9 kW controller power. MATLAB's logged force and power satisfy the published gain equations to numerical precision. The 22.74 m peak-to-peak heave is the published unconstrained linear controller's response; trajectory agreement does not establish a feasible PTO stroke or physical design. |
| Sphere reactive controller with simple direct-drive PTO | Pinned MATLAB Applications `Controls/ReactiveWithPTO`, current Sphere HDF5, and a fresh 200 s ode45 run | The public Python PTO configuration couples the published PI force request to a winding L/R torque state, gear ratio, drivetrain inertia and friction, and generator current/voltage/loss outputs. On 20,001 samples, MATLAB's logged controller, shaft torque decomposition, and electrical identities agree with the source block equations to numerical precision. Its Derivative block logs shaft inertia torque from the preceding output interval's speed difference, matching that finite difference within `3e-10` N m; Python uses instantaneous acceleration in the coupled dynamics. The independent Python trajectory differs by at most 0.076 mm heave and 0.051 mm/s velocity; shaft torque by 0.350 N m and source-signed electrical power by 22.5 W. This validates one pure-heave, zero-heading regular-wave direct-drive connection. Voltage/current limits and other PTO-Sim components remain unsupported. |
| Sphere declutching controller | Pinned MATLAB Applications `Controls/Declutching` and MATLAB-generated Sphere HDF5 | A configurable body-local PTO endpoint uses the published 232,020 N s/m engaged damping and switches off for 0.8 s after velocity reversal. In 40,001 paired samples over 400 s, heave differs by at most 1.99 mm and velocity by 4.89 mm/s. Both runs have 83 disengagement events, each 0.8 s long; event times differ by at most one 0.01 s sample, affecting four force samples. Where both runs are engaged, force differs by at most 1.14 kN and power by 1.06 kW; both report zero force while disengaged. Absorbed energy differs by 10.5 kJ over a 22.4 MJ MATLAB total. A full-sample force maximum is 106 kN because one-step mode differences occur at two transitions, so force parity is gated by event timing and within-mode values. |
| Sphere latching controller | Pinned MATLAB Applications `Controls/Latching` and MATLAB-generated Sphere HDF5 | The public `LatchingControl` applies the published 49,181 N s/m normal damping and 37,308,296 N s/m latch damping for 2.4 s after velocity reversal. Both force laws dissipate energy. In 40,001 paired samples over 400 s, MATLAB and Python each have 83 latch starts; start times differ by at most 0.05 s, and each completed latch lasts 2.41 s in the sampled source implementation. Maximum heave and velocity differences are 0.379 m and 0.573 m/s; the heave extrema differ by at most 0.033 m. Absorbed energy differs by 2.43 MJ against MATLAB's 299.98 MJ (0.81%). This switched high-gain case has 553 samples in different modes, so pointwise controller force is not a useful single parity metric; the paired gate checks the exact force and power laws, event timing, motion, envelope, and integrated energy. The published linear model does not establish a feasible physical stroke or rigid latch. |
| Configured Sphere spring, damper, and PTO attachment | Derived from the pinned Sphere passive case with 50,000 N/m PTO stiffness, 100,000 N s/m native PTO damping, and the PTO moved to x = 1 m | A fresh MATLAB run and the Python `WEC` builder agree within 0.051 mm position and PTO stroke, 0.054 mm/s velocity and PTO speed, 52 N combined force, and 27 W combined power. The x offset has no effect on heave-only motion, so this validates the force settings but not attachment-location dynamics. |
| OSWEC PM equal-energy waves, directional excitation, pitch, and PTO | Current MATLAB OSWEC example and current OSWEC HDF5 | Python recreates all 500 PM equal-energy bins from the HDF5 range, then synthesizes wave elevation and six-component excitation using the saved MATLAB phase matrix. Local component differences are below `6e-15`, elevation below `2e-13` m, and excitation below `1e-7` N. The public two-body `WEC.fixed_hinge` call passed a [fresh pinned R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37795546236) with source phases: maximum flap differences over 4,001 samples were 10.0 mm surge, 1.66 mm heave, 0.00201 rad pitch, 0.00297 rad/s pitch speed, and 35.6 N m PTO torque; the fixed hydrodynamic base agreed exactly. A Python seed produces a separate reproducible realization. This model covers pitch about a fixed hinge, not a general six-DOF device. |
| OSWEC PTO-Sim hydraulic adjustable rod and fixed crank | Pinned MATLAB Applications `PTO-Sim/OSWEC/OSWEC_Hydraulic_PTO` and `OSWEC_Hydraulic_Crank_PTO`, BEMIO OSWEC HDF5, and [fresh 400 s R2025b source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37831809681) | Python couples the published hinged flap and 30 s radiation IRF to either rotary-to-linear linkage, compressible cylinder, rectifying valve, gas accumulators, motor, generator, and PI load. The pinned crank laws map angle to stroke and cylinder force to torque with opposite signs; source-driven cylinder pressure steps, valve flows, accumulator pressures, motor and generator steps agree to numerical precision. Across both independently integrated 40,001-sample trajectories, maximum pitch/speed errors are 0.00136 rad and 0.000722 rad/s, flap-center surge/heave errors 6.55/2.75 mm, crank torque 57.9 kN m against a 2.10 MN m source peak, cylinder force 20.0 kN, chamber pressure 0.390 MPa, accumulator pressure 5.24 kPa, shaft speed 0.0595 rpm, current 0.167 A, and voltage 0.175 V. Valve-port flow differs by up to 0.0540 m³/s at reversals; accumulated port-volume error stays below 0.00129 m³. Explicit paired gates cover body motion, torque, hydraulic states, valve flow, and electrical outputs for each layout. The adjustable-rod source commands negative load resistance after 46 s for 77.0% of its nonzero-current samples; this is a source-case limitation. General hydraulic networks and arbitrary attachment geometry remain unpaired. |
| OSWEC `Multiple_Wave_Spectra` with two independent PM seas | Pinned MATLAB Applications `Multiple_Wave_Spectra` input and BEMIO-generated OSWEC HDF5 | The public `run_case` interface sums separate 2 m, 0° and 1 m, 90° seas with saved independent phase realizations over 100 s. Combined elevation agrees within `1e-11` m and both bodies' six excitation-force components within `1e-6` N or N m. The fixed hydrodynamic base remains stationary. With the ordinary implicit added-mass default, maximum flap differences are 106 mm surge, 10.8 mm heave, and 0.0213 rad pitch; PTO torque differs by at most 7.25 kN m. The explicit source-delay comparison reduces these to 10.3 mm, 1.20 mm, 0.00207 rad, and 818 N m. Both dynamics paths have separate paired gates. The source enables passive-yaw preprocessing, but its hinge has zero yaw motion; the source's spline frequency interpolation is selected explicitly for this case. Moving-yaw dynamics and fixed-base reactions remain unsupported. |
| OSWEC `Full_Directional_Waves` imported spectrum | Pinned MATLAB Applications `Full_Directional_Waves` input and BEMIO-generated OSWEC HDF5 | The imported 47-frequency, 180-heading spectrum and realized phases reproduce all 8,001 wave-elevation samples within `9e-14` m. The source force block omits heading-bin width although its wave-elevation block includes it: logged excitation is about `1/sqrt(2° in radians) = 5.352` times the physically integrated force. Python keeps heading integration as its default. Explicit `matlab_omitted` quadrature plus `spline_frequency` interpolation matches all six logged excitation components for both bodies within `4e-7` N or N m. With that source forcing and ordinary implicit added mass, the 400 s flap differs by at most 12.8 mm surge, 3.52 mm heave, 0.00258 rad pitch, and 56.6 N m PTO torque. The separate opt-in Simulink-delay comparison reduces these to 1.53 mm, 0.293 mm, 0.000306 rad, and 5.54 N m. Source-compatible agreement does not establish that the omitted-width force is physically correct. The fixed base remains stationary; its ground reactions are not computed. |
| OSWEC passive yaw off and on | Pinned MATLAB Applications `PassiveYawOFF` and `PassiveYawON` inputs, BEMIO-generated OSWEC HDF5, and fresh R2025b trajectories | The public `WEC` interface defines a yawing flap, fixed hydrodynamic base, and 120,000 N m s/rad rotational PTO. Both published 600 s cases use 2.5 m, 8 s waves at 10°. With passive yaw off, Python and MATLAB flap yaw agree within `4e-14` rad. With passive yaw on, excitation follows the body-relative wave heading and turns the flap toward the 10° incident heading; yaw and angular-speed differences stay below 0.00149 rad and 0.000193 rad/s. PTO torque and absorbed power differ by at most 23.2 N m and 0.128 W, accounting for MATLAB's negative absorbed-power sign. On the saved MATLAB yaw trajectory, all six reconstructed excitation components differ by at most 100 N or N m except yaw moment, which differs by 467 N m. The source holds heading coefficients until yaw changes by 0.01°, while Python interpolates continuously. The fixed base is stationary; full three-dimensional joint mechanics remain unsupported. |
| OSWEC variable-hydrodynamics passive yaw, 2° bank | Pinned MATLAB Applications `Variable_Hydro/Passive_Yaw` with its published regular-wave, 600 s, 0.01 s case and the first published `-40:2:40` direction bank; fresh R2025b [paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37792909630) | The Python `WEC.body(..., passive_yaw=True, yaw_heading_bank=range(-40, 41, 2))` selects the nearest body-relative BEM heading using the source rule. The selected heading matches MATLAB at all 60,001 saved MATLAB states; yaw excitation evaluated on that path differs by at most 867 N m. Before the first heading-event split at 64.57 s, yaw, yaw speed, and PTO torque differ by at most `3.89e-5` rad, `3.32e-5` rad/s, and 3.98 N m. At 64.57 s the paths differ by only `3.7e-5` rad but select different headings for one sample. Discrete events amplify that offset; over the full run the maximum yaw, speed, and PTO torque differences are 0.164 rad, 0.0221 rad/s, and 2.65 kN m. Full-trajectory parity is **unestablished**. This mode covers the published pure-yaw/fixed-base regular-wave bank, not general variable hydrodynamics. |
| OSWEC irregular passive yaw | Pinned MATLAB Applications `PassiveYawRegression` input, BEMIO-generated OSWEC HDF5, and fresh R2025b trajectory | The public `PMWave` replays the published 500-frequency phase realization over 250 s with 40 s radiation memory. Frequency bins, spectral amplitudes, and widths agree within `6e-15`; elevation agrees within `1e-13` m. The public `passive_yaw_threshold=1` selects the published coefficient hold and nearby BEM-heading snap. On MATLAB's saved yaw path, that opt-in force law matches all six logged excitation components within `1e-7` N or N m. The source's yaw radiation convolution reconstructs from its saved speed and pinned HDF5 within `5.4e-9` N m; added mass, hydrodynamic force balance, and rigid inertia balance close within `1.4e-8` N m. In a native Python run a heading update occurs one 0.01 s sample early at 53.45 s: MATLAB's relative heading is 10.99980° and Python's is 11.00004°, on opposite sides of the 11° update boundary. The yaw difference there is `4.3e-6` rad, and replaying MATLAB's logged force produces nearly the same pre-event drift (`4.29e-6` rad). Later event differences accumulate to 0.203 rad maximum yaw and 0.0277 rad/s speed differences over 250 s. With continuous heading interpolation, Python differs from the published trajectory by up to 0.452 rad yaw and 0.101 rad/s yaw speed. Replaying MATLAB's logged six-component force through the Python dynamics reduces these maximum differences to 0.00180 rad and 0.0000935 rad/s. The source force laws and balances are paired, but published-case trajectory parity with native excitation remains **unestablished** because the threshold amplifies small integration differences. |
| OSWEC irregular passive yaw, continuous-heading control | Same pinned MATLAB input and phase realization, with only `body(1:2).yaw.threshold` changed from 1° to 0° in a temporary application copy | MATLAB and Python then agree over all 25,001 time samples: maximum differences are `1.1e-7` N or N m across the six excitation components on the MATLAB yaw path, `3.16e-5` rad flap/PTO angle, `7.30e-6` rad/s yaw/PTO speed, `0.875` N m PTO torque, and `0.0123` W source-signed power. Wave elevation differs by less than `1e-13` m and the fixed base remains stationary. This control validates the continuous-heading Python dynamics; it does not erase the published 1° threshold difference. |
| OSWEC fixed nonhydrodynamic base | Pinned MATLAB Applications `Nonhydro_Body` case and its BEMIO-generated OSWEC HDF5 | The public `WEC.fixed_hinge(flap, base=wec.fixed_body(...))` call passed a [fresh pinned R2025b paired run](https://github.com/cmudrc/wec-sim-python/actions/runs/37795143695) and reports the stationary base and nonlinear flap motion about the PTO hinge. Across 4,001 samples over 400 s, maximum flap position differences are 20.8 mm surge, 7.9 mm heave, and 0.00444 rad pitch; velocity differences are 15.6 mm/s surge, 6.8 mm/s heave, and 0.00328 rad/s pitch. All six excitation-force components differ by less than 1 N, the base position and velocity agree exactly, and zero PTO torque differs only by MATLAB numerical noise below `3.3e-7` N m. The base's ground-constraint reaction forces are not calculated. |
| Ellipsoid instantaneous nonlinear hydro, `ode4/Regular` | Pinned MATLAB Applications `Nonlinear_Hydro` and BEMIO-generated ellipsoid HDF5/STL | The Python builder combines mesh buoyancy, instantaneous Froude–Krylov correction, quadratic heave drag, BEM linear excitation/radiation, and a configured PTO. On all 3,001 saved MATLAB states, mesh equilibrium mass agrees within `2e-9` kg, buoyancy within `6e-9` N, drag within `7e-11` N, and total heave excitation within `1e-6` N. The independent 150 s Python trajectory differs by at most 5.16 mm heave, 5.86 mm/s velocity, 5.16 mm PTO stroke, 7.03 kN PTO force, and 7.21 kW absorbed power. MATLAB splits the BEM added mass between Simscape mass and an applied force; Python uses the combined effective mass. This case validates the single-body, pure-heave, zero-direction constant-radiation mode. |
| Ellipsoid instantaneous nonlinear hydro, `ode4/RegularCIC` | Pinned MATLAB Applications `Nonlinear_Hydro` RegularCIC input and BEMIO-generated ellipsoid HDF5/STL | The same Python mesh model uses a 60 s radiation convolution. On all 3,001 saved MATLAB states, mesh buoyancy and quadratic drag agree within `6e-9` and `7e-11` N, linear plus nonlinear heave excitation within `4e-7` N, and radiation convolution within `3e-10` N. The independent 150 s Python trajectory differs by at most 6.39 mm heave, 7.41 mm/s velocity, 6.39 mm PTO stroke, 8.89 kN PTO force, and 8.98 kW absorbed power. The memory-step fixed-point iteration was allowed more iterations to converge near the moving waterline; no force coefficient was tuned. |
| Ellipsoid `ode45/Regular` and `ode45/RegularCIC` diagnostics | Pinned MATLAB Applications `Nonlinear_Hydro` ode45 inputs and BEMIO-generated ellipsoid HDF5/STL | The physical inputs match the two ode4 cases, but MATLAB ode45 applies nonlinear buoyancy from the preceding 0.05 s output sample: its restoring log matches the mesh law on the prior body state within `6e-9` N, while the current-state mismatch exceeds 54 kN. Its logged total force and adjusted Simscape mass reconstruct the logged acceleration within `2e-9` N, confirming that this is an applied source force. Drag and wave excitation remain current-state forces; the RegularCIC radiation log differs from convolution of output-sampled velocity by at most 266 N. MATLAB ode4 versus ode45 differs by up to 15.9 mm heave and 35.4 kN PTO force with identical physical settings. Python's instantaneous-force trajectories differ from MATLAB ode45 by at most 19.94 mm heave, 30.05 mm/s velocity, 36.05 kN PTO force, and 28.78 kW absorbed power. These are solver-envelope diagnostics, **not established ode45 numerical parity**; Python does not add a one-sample buoyancy delay to reproduce the source solver artifact. |
| Sphere variable draft and mass in heave | Pinned MATLAB Applications `Variable_Hydro/Variable_Mass`, nine BEMIO-generated draft HDF5 files, and its full 900 s ode4 run | The Python `WEC.variable_body` selects the published state every 100 s and changes equilibrium draft, rigid mass, added mass, radiation damping, restoring stiffness, and excitation together. On MATLAB's 9,001 saved states, state-selected excitation differs by at most `3e-8` N, radiation damping by `3e-10` N, restoring plus weight by `8e-9` N, and assembled hydrodynamic force by `9e-9` N. The independent Python trajectory differs by at most 2.23 mm heave, 4.14 mm/s velocity, and 827 N PTO force over all nine states. The source reaches −6.63 m by 900 s; agreement with its linear coefficients does not establish physical accuracy at that excursion. Other variable-hydrodynamic motions remain unsupported. |
| Barge generalized body modes | Pinned MATLAB Applications `Generalized_Body_Modes`, its BEMIO-generated ten-DOF barge HDF5, and a fresh 400 s R2025b run including the model's separate `Flex_out` signal | `BodyClass` pairs flexible effective mass, stiffness, damping, hydrostatic, added-mass, radiation-damping, and excitation matrices; the largest absolute matrix difference is `1.91e-6` in applied added mass. The source flexible acceleration equation and six signed force components reconstruct from those matrices, with restoring and radiation differences below `1e-6` N, and sampled added-mass feedback below `0.05` N. The public `WEC.floating_gbm` interface runs a coupled three-rigid/four-flexible-coordinate implicit-mass model in the published 2 m, 8 s regular wave. Over all 8,001 samples, maximum rigid position and velocity differences are 1.06 mm and 0.838 mm/s; flexible displacement, velocity, and acceleration differ by at most 0.043 mm, 0.034 mm/s, and `1.26e-4` m/s². This validates one single-body floating 3-DOF joint with four modes, zero-heading regular waves, and no PTO or mooring; other GBM layouts and wave/radiation options remain unsupported. |
| OWC `OrificeModel` wave, force, and early coupled motion | Pinned MATLAB Applications `OWC/OrificeModel`, BEMIO-cleaned seven-DOF HDF5, and a fresh 130 s R2025b `ode23t` [source run](https://github.com/cmudrc/wec-sim-python/actions/runs/37799665808) | `wecsim.OrificePTO` reproduces the embedded Simulink orifice block on its 418,828 logged internal solver evaluations. Maximum force, pressure, and source-power differences are `4.8e-10` N, `2.9e-12` kPa, and `4.2e-13` kW; the Mach flag is identical. Interpolating source flow to the independent 26,001-sample flexible-mode output gives a piston-speed difference below `3.6e-13` m/s. The 500 PM bins agree within `5.2e-15`, wave elevation within `2.8e-14` m, and all six rigid plus one flexible excitation channels within `1.5e-10` N. Flexible mass/stiffness/damping, added mass, and hydrostatic matrices agree within `2.3e-13` in their units. With the [rigid-force diagnostic run](https://github.com/cmudrc/wec-sim-python/actions/runs/37804483371), radiation reconstructed from all seven MATLAB velocities differs by at most `0.0074` N in the flexible channel and `0.22` N across rigid channels (peak rigid magnitude `1,182` N). The rigid restoring and linear-damping forces reproduce the source logs within `8.6e-11` N and `1.1e-12` N, and the rigid hydrodynamic force sum closes within `1.7e-10` N. MATLAB's flexible hydrodynamic force sum closes within `1e-12` N; its effective-mass acceleration including the orifice force closes within `3.5e-8` N. The source rigid surge, heave, and pitch balances close within `1.2e-10` in their force units when using its adjusted translational mass and logged piston reaction. `WEC.floating_gbm(..., orifice=...)` now integrates the seven-DOF PM case with memory radiation, quadratic drag, and the coupled piston force. Against the saved source phases, maximum surge, heave, unwrapped-pitch, and flexible-position differences through 6 s are `0.0113` m, `0.0244` m, `0.0098` rad, and `0.0070` m. At 130 s the heave and flex differences reach `1.94` m and `0.280` m; **full trajectory parity is not established**. Python reports SI watts while the source's numeric power is kW. The source first exceeds its Mach threshold at 10.13 s; it is above 0.3 for 61.4% of simulated time and peaks at Mach 1.11. Matching its incompressible law after that point does not establish physical accuracy. |
| Case-driven dynamics runner | Current MATLAB RM3 and OSWEC examples plus Sphere and RM3 Applications cases | One generalized-coordinate engine assembles rigid inertia, hydrodynamic added mass and radiation, hydrostatic restoring, excitation, and linear PTO forces. The heave, fixed-hinge, and floating-joint layouts run the paired cases above, including imported elevation and a joint surge spring for RM3 MooringMatrix. A fourth `linear_subspace` layout maps independent coordinates into arbitrary bodies; its RM3 two-body heave, Sphere free decay, configured Sphere PTO, and two instantaneous nonlinear-hydro cases have paired MATLAB checks. Arbitrary Simscape layouts, general moorings, and other nonlinear-hydro/application cases remain unsupported. |

The pinned OWC Simulink model has a force-path difference:
its flexible `GBM` state-space block integrates inverse effective mass times
the hydrodynamic force, while the piston force joins a separate reported
acceleration path. Over the first 6 s, a finite difference of the saved
flexible velocity agrees with hydrodynamic force divided by effective mass to
`0.00196` m/s² RMS; it differs from the reported acceleration that includes
the piston by `0.257` m/s² RMS. At 2 s these three values are `-0.00317`,
`-0.00299`, and `-0.15217` m/s², respectively. The rigid-heave added-mass path
still receives the reported flexible acceleration. The opt-in
`orifice_force_path="published_owc"` uses those published signal routes and
reduces the 130 s heave and flexible
position differences from `1.94` m and `0.280` m to `0.0366` m and `0.00533` m.
This identifies the primary cause of their divergence. Python keeps the piston
reaction in the flexible state equation by default. The source-compatible
setting is for trajectory comparison; it does not restore the omitted piston
force or establish physical validity above the source's Mach threshold. Surge
and unwrapped pitch still differ by up to `0.174` m and `0.0107` rad, so exact
full-trajectory numerical parity remains unestablished.

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
The diagnostic fixture and a freshly BEMIO-generated HDF5 from the pinned
RM3 Applications input have identical 260 frequency samples and source BEM
damping. Their common-surge fitted transfer curves differ by at most
`6.9e-9` N s/m and projected active-mode minimum damping by at most
`3.0e-9` N s/m. The raw fitted A/B/C arrays differ because they encode
equivalent state-space realizations; the paired gate compares the transfer
responses. Thus the negative low-frequency fit is present in the current
MATLAB case input as well as the historical fixture.
The diagnostic also projects the BEM and fitted radiation matrices onto
the RM3 joint's four moving coordinates: common surge, float heave, spar
heave, and shared pitch. Pitch is expressed as travel at a 20 m lever for
an energy-conjugate matrix with consistent units. The fitted matrix's
least-damped direction is −14.3 kN s/m at zero frequency with cross-body
radiation off and −17.6 kN s/m with it on. The source BEM matrix has small
negative eigenvalues too (minimum −0.422 kN s/m across its sampled
frequencies), so these data do not establish global BEM passivity. At the
sampled frequency nearest the 8 s wave, the coupled source minimum is
−0.030 kN s/m while the fit is −5.87 kN s/m. These are properties of the
linearized radiation data in the active joint subspace, not proof that every
negative mode is excited in Cases 5 and 6.
Diagnostic PR #6 was closed without merging its state-space solver. The
convolution results for Cases 3 and 4 were identical, sample for sample,
between the original branch point and that diagnostic branch.
The pitched-slider map used by the current default solver is also checked
independently of MATLAB output: finite differences of body position verify
its velocity Jacobian and acceleration bias, and an unforced two-body joint
with symmetric added mass conserves kinetic energy within `1e-8` relative
over 10 s. These checks guard the default joint mechanics; they do not
establish passivity of the published fitted radiation model.
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

The [latest full reference-model sweep](https://github.com/cmudrc/wec-sim-python/actions/runs/37840433084)
completed 34 MATLAB/Python jobs successfully, with one optional comparison
skipped. That sweep covers the cases selected by its matrix, not all 63
inventory entries. The [separate current-profile run](https://github.com/cmudrc/wec-sim-python/actions/runs/37820035413)
covers the uniform, power-law, and linear fixed-Morison comparisons added
after the sweep's code revision.

### Named gaps in the published case inventory

The table identifies cases with incomplete independent Python motion/output
coverage at the pinned source revisions. The detailed paired limits for covered
cases are above. A passing upstream MATLAB application test only establishes
that WEC-Sim ran; it does not establish Python parity.

| Published case(s) | Evidence and remaining work |
| --- | --- |
| RM3 `B2B_Case5` and `B2B_Case6`; state-space setting in `Radiation_Force_Options` | MATLAB trajectories and fitted-transfer diagnostics exist. The pinned fit has negative low-frequency damping in the active joint coordinates. No physically validated Python state-space trajectory pair exists. |
| OSWEC `PassiveYawRegression` and `Variable_Hydro/Passive_Yaw` | Source force laws and early motion are paired. Discrete heading updates amplify small integration offsets, so published full-trajectory parity is unestablished. |
| RM3 `Mooring/MoorDyn` and `Paraview_Visualization/RM3_MoorDyn_Viz` | The pinned upstream CI skips MoorDyn; forcing it on GitHub exposed a native-library mismatch. There is no paired Python mooring trajectory. |
| OSWEC `Desalination` | The pinned 300 s Simscape Fluids application runs under R2025b; [expanded source baseline](https://github.com/cmudrc/wec-sim-python/actions/runs/37838787789) records flap, PTO, pressure, flow, and mechanical-power traces. Its 250-bin incident sea, six flap excitation components, and body-local rod motion pair with Python over all 30,001 samples; maximum excitation error is below `1e-6` N or N m. The source Morison drag law matches 605 sampled flap states within `1e-6` N or N m after accounting for the source's logged-force sign. The source PTO actuation force exactly follows the measured cylinder force with a one-step delay over the full trajectory. Given MATLAB's inlet pressure, the Python osmotic valve and resistance reproduce all 30,001 permeate-flow values within `2.1e-12` m³/s. Given the source accumulator flow and startup pressure, the Python gas-volume and hard-stop law reconstructs all 30,001 pressure values within `1e-3` Pa using the source's backward Euler flow step. The initial inferred liquid volume is negative because the source starts its network far below the 3 MPa precharge while allowing finite hard-stop penetration. Given the source chamber pressures and rod speed, the Python ideal cylinder reproduces all 30,001 rod-force samples within `1e-7` N and both port-flow traces within `1e-12` m³/s; the measured source force has the opposite sign. The published hydraulic junction balance closes within `2e-15` m³/s, including the relief branch. The high-pressure network now advances from saved rod speed alone, without measured source pressure or branch flows. Over all 30,001 samples it gates pressure within `500` Pa (source pressure reaches 5.6 MPa), permeate, brine, and recovered flow within `1e-5` m³/s each, and relief/accumulator flow within `3e-3` m³/s each. Using the source feed flow instead, the network reconstructs the first 6001 pressure states within `0.01` Pa. This also provides a stronger prescribed-rod-motion component gate than the coupled check below. A four-valve cylinder model using the published passive-orifice settings predicts chamber pressures and rod force from Python high pressure and the saved rod speed. The pinned legacy solver alternates raw chamber pressure at the 0.01 s sample rate while valve flow remains smooth; over 20% of valve 1 samples log flow against the pressure drop. Pointwise chamber pressure and force are therefore not claimed as paired. After a 0.1 s mean, the independent rod-force error is below 100 kN RMS; net rod work differs by less than 3% over 300 s. Raw rod-force disagreement remains above 1 MN RMS. The independent 300 s Python runner now couples incident sea, flap motion, radiation memory, source-convention Morison drag, four-valve chamber force, and high-pressure dynamics. With the pinned HDF5 and saved MATLAB random phases, but no source forces or state trajectories as inputs, its paired gates are maximum flap-pitch error below `0.02` rad, pitch-speed error below `0.01` rad/s, body-center error below `0.1` m, rod-speed error below `0.02` m/s, high-pressure error below `40` kPa, rod-stroke error below `0.05` m, 0.1 s mean rod-force error below `150` kN RMS, and net rod-work error below `2%` over 300 s. The source-specific Morison convention and one-step PTO actuation delay are confined to this case runner. The raw chamber pressure/force artifact remains unpaired. The source also permits chamber gauge pressure below −15 MPa, so this numerical comparison is not evidence of physical cavitation behavior. The two published OSWEC PTO-Sim hydraulic crank applications are covered above. |
| OSWEC `Paraview_Visualization/OSWEC_NonLinear_Viz`; RM3 and OSWEC `Wave_Markers` | Their visualization or marker outputs have no direct Python comparison. Regular-wave body dynamics are paired elsewhere but do not verify these outputs. |
| Other/dynamic: `Load_Mitigating_Controls/ControlTests`, `MOST`, `OWC/FloatingOWC`, `WECCCOMP/WECCCOMP`, `WECCCOMP/WECCCOMP_Fault_Implementation`, and `WECCCOMP/WECCCOMP_Nonlinear_Model_Predictive` | Inventoried inputs without a published-case Python motion/output gate. They require case-specific mechanics or controls. The preceding WaveBot `CalcImpedance` system-identification stage is paired separately. |
| Other/dynamic: both `Nonlinear_Hydro/ode45` cases and `OWC/OrificeModel` | Source force components and bounded motion intervals are compared above, but full-trajectory parity remains unestablished because of source solver force timing or large late motion differences. |

The published RM3 `End_Stops` trajectory has a separate refined-step gate;
the pinned 0.1 s output is not treated as a converged motion reference.

The `MATLAB reference model baselines` workflow runs the two canonical core
examples, all five Sphere free-decay cases, RM3 body-to-body Cases 1–4,
the RM3 MooringMatrix case when selected by workflow dispatch,
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
remains available for saved cases; the Python builder covers the mapped
`linear_subspace`, one-body `floating_gbm`, two-body RM3-style
`floating_joint`, and one- or two-body OSWEC `fixed_hinge` layouts. The
floating joint uses a relative-heave PTO; the fixed hinge accepts a torsional
PTO at the published world-axis locations. Neither accepts arbitrary
attachment geometry.
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
Separate direct runs of the pinned [published irregular passive-yaw input](https://github.com/cmudrc/wec-sim-python/actions/runs/37718464539)
and its [0° heading-update control](https://github.com/cmudrc/wec-sim-python/actions/runs/37720514831)
completed and generated the paired data summarized above. Those runs do not
change the upstream regression assertions.

## Known differences and next reference case

Current MATLAB WEC-Sim changed its PM/JS spectra, seeded phase generator,
object properties, and some wave inputs since the original Python port.
The general `WaveClass` now has direct paired Traditional and EqualEnergy
PM/JS checks using replayed MATLAB phases; the focused OSWEC PM and other
application-specific irregular-wave paths above have separate paired gates.
The historical fixtures alone do not establish parity for remaining wave modes
or independent seeded realizations across MATLAB and NumPy.
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
