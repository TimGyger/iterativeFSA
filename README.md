# Iterative Methods for Full-Scale Gaussian Process Approximations

This repository provides the R code for the simulation studies and real-world applications presented in the paper *"Iterative Methods for Full-Scale Gaussian Process Approximations for Large Spatial Data."*  

The iterative methods for full-scale approximations are implemented in the **GPBoost** package, available here: [https://github.com/fabsig/GPBoost](https://github.com/fabsig/GPBoost).

## Setup

This repository uses Git submodules for external dependencies.

For cloning, run:
```
git clone --recurse-submodules https://github.com/TimGyger/iterativeFSA.git
```

If you already cloned the repository, initialize submodules with:
```
git submodule update --init --recursive
```

## Repository Structure

- **`Data`**  
  Provides scripts for generating data for both the simulation studies and real-world experiments.

- **`Simulation Studies`**  
  Contains scripts to run the simulation experiments. Refer to the `README.md` file in the `simulation_studies` folder for a detailed description of the various simulations and how to execute them.

- **`Real World Application`**  
  Includes scripts for the real-world application on daytime land surface temperature data.
