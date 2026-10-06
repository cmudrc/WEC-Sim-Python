# WEC-Sim-Python

> **cmudrc fork status:** This is an active parity effort, not yet a complete
> wave energy converter simulator. A focused linear heave free-decay solver
> is available; the original general runner still has no dynamics stage.
> See [PARITY.md](PARITY.md) for verified behavior, current
> MATLAB reference revision, and the remaining work.

To run the production-code parity checks with Python 3.12:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt pytest
.venv/bin/python -m compileall -q source/objects
.venv/bin/python -m pytest -q tests/test_wave_parity.py tests/test_body_io.py
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

response = solve_heave_free_decay("sphere.h5", initial_displacement=1.0)
# response.time, response.position, response.velocity, response.force_total
```

This solver assumes one heave-only body, zero incident waves, and no PTO,
mooring, or nonlinear force. The original README projected completion in
August 2022; that date is no longer applicable. See [PARITY.md](PARITY.md)
for the tested scope and next reference case.

## Design direction

The Python implementation will solve device dynamics without Simulink. Its
preprocessing, equations, and numerical outputs must be compared to explicit
MATLAB WEC-Sim revisions and reference cases before a feature is called
supported. ParaView file output and BEMIO conversion are future work.
