# Thermal Protection Simulation (1D Multi-Layer Heat Conduction)

## Overview

This program simulates 1D transient heat conduction across a stack of thermal protection materials. The robot's body is modeled along the vertical axis (z-axis), and each layer's thermal properties (e.g., conductivity, specific heat) are considered. The simulation calculates the temperature profile over time and outputs a full dataset for further analysis.

## Features

- Solves the 1D heat equation using the **Backward-Time Central-Space (BTCS)** method (implicit).
- Models varying layer thicknesses along the robot's body (`z` from 0 to 2.5 m).
- Includes temperature-dependent boundary conditions.
- Outputs:
  - Minimum required insulation thickness at each `z`.
  - Full temperature profiles over time and space into CSV.

## File Outputs

- `Thickness_<duration>_hr.csv`: Required insulation thickness vs. height.
- `full_temp_profile.csv`: Time-resolved temperatures across all materials and `z`.

## Material Stack

Each node is composed of four layers:
1. **Thermal Protection** (externally exposed, thickness computed)
2. **Carbon Fiber**
3. **Glue**
4. **Steel**

Each has configurable thermal properties and resolution (number of grid points per layer).

## Dependencies

- C++17 or newer (for `std::string`, `vector`, `cmath`)
- No third-party libraries

---

## How It Works

### Step 1: Initial Thickness Calculation

```cpp
vector<double> calculateRequiredInsulationThicknessMultiLayer(...);
```

- Runs a binary search for the minimum thickness needed to keep the carbon fiber below its glass transition temperature.
- Solves for steady-state temperature after `missionDuration`.

### Step 2: Time Evolution

```cpp
void MakeTimeSolution(...);
```

- Simulates the full stack for all `z` values and timesteps.
- Outputs the layer-wise temperature for each time and height into `full_temp_profile.csv`.

### Step 3: Multilayer Solver

```cpp
vector<vector<double>> solveMultiLayer(...);
```

- Builds the global temperature grid.
- Handles interface nodes using harmonic means of thermal diffusivities.
- Solves with the **Thomas algorithm** for each timestep.

---

## Data Output Format (`full_temp_profile.csv`)

Each row is formatted as:

```csv
t=<time>,z=<z>,<LayerName1>,<Thickness>,<dx>,T0,...,TN,<LayerName2>,...
```

### Example:

```csv
t=0,z=0.1,THERMAL,18.0,0.002,...,CARBON,1.4,0.001,...,GLUE,...
```

---

## Constants and Functions

- `initialTemp(z)`: External boundary temperature (logarithmic profile).
- `carbonThickness(z)`: Sinusoidal variation
- `glueThickness(z)`: Logarithmic variation
- `steelThickness(z)`: Sawtooth variation
- `MaterialProperties`: Holds `k`, `ρ`, `cp`, `alpha`, `name`

---

## Configurable Parameters

You can modify these in `main()`:

```cpp
double dz = 0.1;              // vertical resolution
double dt = 1.0;              // time step
double missionDuration = 3600; // seconds
vector<int> pointsPerLayer = {10, 10, 10, 10}; // resolution per layer
```
---

## 🧪 Testing

Compile with:

```bash
g++ -std=c++17 -o thermal_solver 1D_heat_solver_final.cpp
./thermal_solver
```
