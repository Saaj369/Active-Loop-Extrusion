<div align="center">

# 🧬 Active Loop Extrusion of Chromatin

**MS Thesis · Computational Biophysics**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: C++](https://img.shields.io/badge/Language-C%2B%2B-00599C?logo=c%2B%2B)](https://isocpp.org/)
[![Analysis: Python](https://img.shields.io/badge/Analysis-Python-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Langevin Dynamics](https://img.shields.io/badge/Method-Langevin%20Dynamics-8A2BE2)]()

*Simulating and studying the  Active loop extrusion hypothesis in chromatin organization using Langevin dynamics*

---


<!-- Replace the src below with your actual animation GIF paths once uploaded -->
<!-- <img src="assets/animations/ring_animation.gif" alt="Ring Polymer Loop Extrusion" width="45%"/>
&nbsp;&nbsp;
<img src="assets/animations/linear_animation.gif" alt="Linear Polymer Loop Extrusion" width="45%"/>

*Left: Ring polymer undergoing active loop extrusion &nbsp;|&nbsp; Right: Linear polymer dynamics* -->

---

</div>

## Overview

This repository contains all simulation and analysis code developed for my MS thesis on **Active Loop Extrusion of Chromatin**. The work investigates how molecular motors (such as cohesin and condensin) actively extrude chromatin loops, driving the spatial organization of the genome inside the cell nucleus.

Loop extrusion is modeled as an active, non-equilibrium process on a polymer chain, with the extruder performing directed motion along the chain while anchored at both feet. The thesis explores two polymer topologies:

- **Ring Polymer** — a closed circular chain; serves as the foundational model (overdamped regime only)
- **Linear Polymer (Realistic)** — an open-ended chain mimicking actual chromatin; studied in both overdamped and underdamped regimes

Simulations are written in **C++** for performance, and post-processing, analysis, and visualization are done in **Python**.

---

## 🗂️ Repository Structure

```
Active-Loop-Extrusion/
│
├── Ring Polymer/                    # Closed-chain topology model
│   ├── Overdamped/                  # Overdamped Langevin dynamics
│   │   ├── src/                     # C++ simulation source files
│   │   ├── analysis/                # Python analysis & plotting scripts
│   │   └── readme.txt               # Build & run instructions
│   └── ...
│
├── Linear Polymer (Realistic)/      # Open-chain (realistic chromatin) model
│   ├── Overdamped/                  # Overdamped Langevin dynamics
│   │   ├── src/
│   │   ├── analysis/
│   │   └── readme.txt
│   ├── Underdamped/                 # Underdamped (inertial) Langevin dynamics
│   │   ├── src/
│   │   ├── analysis/
│   │   └── readme.txt
│   └── ...
│
├── assets/
│   ├── animations/                  # GIF/MP4 visualizations of simulations
│   └── thesis.pdf                   # Full MS Thesis document
│
└── README.md
```

> 📄 Each subdirectory contains a `readme.txt` with detailed instructions on how to **build and run** the corresponding simulation.

---

## 🔬 Physical Models

### 1 · Ring Polymer Model

https://github.com/user-attachments/assets/a21c97d1-8cdb-4200-b269-249d957c33fd

The bacterial chromatin is generally arranged in circular structure therefore Ring Polymer was a natural starting point. The ring polymer is modeled as a bead-spring chain forming a closed loop. The extruder binds on a random site of polymer and pulls the polymer in. This causes the ring to collapse and get folded.
- **Regime:** Overdamped Langevin dynamics
- **Integration scheme:** Velocity Verlet
- **Key observables:** Loop size distribution, mean-square displacement, contact maps

### 2 · Linear Polymer (Realistic) Model


https://github.com/user-attachments/assets/03d34c3c-0eb6-4a9c-8903-a0f09f08fab3



The linear model more faithfully represents chromatin as a linear chain with free ends, as found in eukaryotic chromosomes. This model introduces CTCFs as boundary elements. The novelty of this model is that it attempts to create realistic model of a chromatin by placing CTCF at polymer locations identified from the real genomic data. This model was studied in two dynamical regimes:

| Regime | Description |
|---|---|
| **Overdamped** | Friction dominates; inertia neglected. Appropriate for highly viscous nucleoplasm. |
| **Underdamped** | Inertia retained; captures transient dynamics and momentum effects at shorter timescales. |

---

## ⚙️ Simulation Details

| Property | Value |
|---|---|
| Simulation engine | Custom C++ (compiled with `g++`) |
| Integration method | Langevin (Velocity Verlet) |
| Polymer model | Bead-spring (FENE / Harmonic bonds) |
| Extruder dynamics | Active directed stepping + stochastic unbinding |
| Analysis & plotting | Python (`NumPy`, `Matplotlib`, `SciPy`) |

---

## 🚀 Getting Started

Each model directory contains a `readme.txt` with specific build instructions. The general workflow is:

### 1. Clone the repository

```bash
git clone https://github.com/Saaj369/Active-Loop-Extrusion.git
cd Active-Loop-Extrusion
```

### 2. Build the C++ simulation

Navigate to the desired model directory and follow the instructions in `readme.txt`. A typical build looks like:

```bash
cd "Ring Polymer/Overdamped"
g++ -O3 -o simulate src/main.cpp -lm
```

### 3. Run the simulation

```bash
./simulate
```

Output trajectory/data files will be generated in the working directory.

### 4. Analyze and plot

```bash
cd analysis/
python3 plot_results.py
```

> ⚠️ Please refer to the `readme.txt` inside each subdirectory for exact compiler flags, parameter files, and expected output formats.

---

## 🎬 Animations

Visualizations of the polymer dynamics and loop extrusion process are available in the [`assets/animations/`](assets/animations/) directory.

<!-- Add more rows as needed for your animations -->
| Animation | Description |
|---|---|
| `ring_animation.gif` | Loop extrusion on a ring polymer |
| `linear_overdamped.gif` | Linear chain — overdamped regime |
| `linear_underdamped.gif` | Linear chain — underdamped regime |

---

## 📄 Thesis

The full MS thesis document is available here:

**[📥 Download Thesis (PDF)](assets/thesis.pdf)**

> *Title:* **Active Loop Extrusion of Chromatin**
> Contains detailed theoretical background, model derivations, simulation methodology, results, and discussion.

---

## 🧑‍💻 Dependencies

### C++ (Simulation)
- `g++` (C++17 or later recommended)
- Standard library only — no external dependencies

### Python (Analysis & Plotting)
- `numpy`
- `matplotlib`
- `scipy`
- *(optional)* `pandas`, `seaborn` for extended analysis

Install Python dependencies with:

```bash
pip install numpy matplotlib scipy
```

---

## 📬 Contact

Feel free to reach out for questions, discussions, or collaborations.

**Saaj** — *MS in Physics / Computational Biophysics*
[![GitHub](https://img.shields.io/badge/GitHub-Saaj369-181717?logo=github)](https://github.com/Saaj369)

---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

<div align="center">

*If this work is useful to you, consider leaving a ⭐ on the repository!*

</div>
