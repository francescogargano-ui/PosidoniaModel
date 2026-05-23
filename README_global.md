# PosidoniaModel

This repository contains the numerical codes and post-processing tools used for the computational study of spatial self-organisation and resilience in *Posidonia oceanica* meadows.

The repository is organized around three main components:

1. a 2D FEniCS simulation code for the reduced spatial model;
2. a Julia/BifurcationKit code for the 1D bifurcation analysis;
3. documentation for accessing the MATLAB `.fig` source files and extracting figure data.

## Repository contents

| File | Description |
|---|---|
| `fenics_code.py` | Clean Python/FEniCS implementation of the 2D reduced *Posidonia oceanica* model. It performs the nondimensionalisation, builds the computational mesh, sets the initial condition, solves the PDE system, and saves the numerical output. |
| `fenics_code.md` | Short documentation for running `fenics_code.py` and visualizing the saved density field in MATLAB. |
| `posidonia_bifurcation.jl` | Clean Julia code for the 1D bifurcation analysis. It defines the model, builds the finite-difference operators, continues the homogeneous branch, detects the first bifurcation, and continues one secondary branch. |
| `Bifurcation_Code.md` | Short documentation for the Julia bifurcation code and for the JSON export of computed branches. |
| `plot_branch_json.m` | MATLAB function used to plot bifurcation branches saved in JSON format by the Julia code. |
| `Data_source_figures.md` | Instructions for accessing the MATLAB `.fig` source files and extracting numerical data from the figures. |
| `LICENSE` | MIT license for the code distributed in this repository. |

## 2D FEniCS simulations

The file

```text
fenics_code.py
```

runs the 2D finite-element simulations of the reduced spatial model.

The code includes:

- nondimensionalisation of the model parameters;
- construction of the 2D computational domain;
- initialization with localized vegetation patches;
- mixed finite-element formulation for the shoot density and its Laplacian;
- time integration with fixed time step;
- XDMF/HDF5 output for visualization and post-processing.

Run with:

```bash
python fenics_code.py
```

The detailed notes are provided in:

```text
fenics_code.md
```

## 1D bifurcation analysis

The file

```text
posidonia_bifurcation.jl
```

performs the 1D bifurcation analysis using Julia and BifurcationKit.

The code:

- defines the 1D reduced model;
- constructs finite-difference operators with no-flux boundary conditions;
- continues the homogeneous branch with respect to the mortality parameter;
- detects the first bifurcation point;
- continues one secondary branch emerging from that bifurcation;
- saves the computed branches in JSON format.

Run with:

```bash
julia posidonia_bifurcation.jl
```

The detailed notes are provided in:

```text
Bifurcation_Code.md
```

## Plotting bifurcation branches

The Julia code saves bifurcation branches as JSON files. These can be plotted in MATLAB using:

```text
plot_branch_json.m
```

Example:

```matlab
figure; hold on;

plot_branch_json('data/br_primary.json', 'meanu');
plot_branch_json('data/br_from_first_bp.json', 'meanu');

xlabel('\omega_{d0}');
ylabel('mean(n)');
title('Bifurcation diagram');

hold off;
```

Other saved quantities, such as `maxn` or `amp`, can also be plotted by changing the second argument.

## Figure source files

The MATLAB `.fig` files used to generate the manuscript figures are not necessarily stored directly in the repository, because some of them may be too large for standard GitHub upload.

Instructions for accessing, opening, and extracting data from the figure source files are provided in:

```text
Data_source_figures.md
```



```
