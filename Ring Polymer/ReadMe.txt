Ring Polymer Active Loop Extrusion Simulation
=============================================

This repository contains C++ code for simulating **Active Loop Extrusion** on a **ring polymer**, combining elements of molecular dynamics (MD) and stochastic SMC (Structural Maintenance of Chromosomes) behavior.

The code was developed as part of a research project modeling chromatin folding through loop extrusion mechanisms.

------------------------------------------------------------
Structure of the Code
------------------------------------------------------------

Main simulation:
----------------
- `SMC_Ring.cpp` : Entry point of the simulation. Integrates polymer dynamics and SMC actions on a ring polymer.
- `mkSMC`        : Makefile to compile the simulation.

SMC logic:
----------
- `SMC.cpp`, `SMC.h`               : Handles SMC-related logic — binding, unbinding, loop extrusion.
- `SMC_helpers.cpp`, `SMC_helpers.h` : Utility functions supporting SMC actions.

Molecular Dynamics (MD) components:
-----------------------------------
- `helpers.cpp`, `helpers.h`       : Core helper functions for MD simulation (forces, integration).
- `helpers2.cpp`, `helpers2.h`     : Extended helpers for additional features or performance tuning.
- `parameters.h`                   : Global simulation parameters (e.g., temperature, time step, number of monomers).

------------------------------------------------------------
How to Build
------------------------------------------------------------

1. Make sure you have a C++ compiler (g++ recommended).
2. Run the Makefile: make -f mkSMC


3. This will compile and generate an executable (e.g., `SMC_Ring`).

------------------------------------------------------------
Usage & Output
------------------------------------------------------------

- Run the executable to perform a loop extrusion simulation.
- Output data (positions, SMC states, etc.) will be saved as `.dat` or `.txt` files.
- Parameters can be modified in `parameters.h`.

------------------------------------------------------------
Citation & Context
------------------------------------------------------------

This code is part of the thesis project:

**"Chromosomal Compaction and Domain Formation by Active Loop Extrusion"**

If you use or adapt this code, please cite or acknowledge appropriately.

------------------------------------------------------------
Author
------------------------------------------------------------

Sajan Daheriya (Saaj)  
Physics Major, IISER Bhopal

