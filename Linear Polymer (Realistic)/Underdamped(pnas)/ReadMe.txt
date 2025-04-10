Linear Polymer Simulation with Underdamped Langevin Dynamics
=============================================================

This code simulates **active loop extrusion on a linear polymer**, using **underdamped Langevin dynamics**. It combines polymer physics with active SMC dynamics, mimicking chromosomal organization processes such as loop extrusion by cohesin-like motors.

The simulation uses a coarse-grained polymer model where each monomer represents a segment of DNA.

------------------------------------------------------------
Structure of the Code
------------------------------------------------------------

Main file:
----------
- `DNA.cpp` : All-in-one simulation file. Contains dynamics, SMC logic, data output, and initialization.  
              Uses OpenMP for parallel force updates and a fast RNG for stochastic components.

Dependencies:
-------------
- `<p2rng/p2rng.hpp>` : Header-only RNG used in place of standard `rand()` for better randomness.
- Standard C++ headers: `<iostream>`, `<cmath>`, `<random>`, `<fstream>`, etc.
- `parameters.h` : Stores simulation parameters (number of monomers, time step, temperature, etc.).
- `monomer.h`    : Defines the monomer structure and attributes (positions, velocities, etc.).

Note: The code is written in a modular way inside one main file (`DNA.cpp`) for simplicity and compactness.

------------------------------------------------------------
How to Build
------------------------------------------------------------

Requires a C++11 (or later) compatible compiler.

1. Compile using g++ with OpenMP enabled:

g++ -O3 -fopenmp -std=c++11 DNA.cpp -o DNA_sim


2. Run the executable:

./DNA_sim


------------------------------------------------------------
Output
------------------------------------------------------------

- The simulation will produce output files for:
- Monomer positions over time
- Active SMC positions / loop states
- Output format can be customized in the code (CSV / tab-separated text, etc.).

------------------------------------------------------------
Project Context
------------------------------------------------------------

This simulation is part of the project:

**"Chromosomal Compaction and Domain Formation by Active Loop Extrusion"**

Linear polymer model complements the ring polymer simulation by modeling boundary effects and open-ended DNA structure.

------------------------------------------------------------
Author
------------------------------------------------------------

Sajan Daheriya (Saaj)  
Physics Major, IISER B

