# Repeatability Intructions

## Requirements

- MATLAB >R2017
- SReachTools [https://unm-hscl.github.io/SReachTools/](https://unm-hscl.github.io/SReachTools/)
- Model Parametric Toolbox [MPT3.0](https://www.mpt3.org/)
- CVX [http://cvxr.com/cvx/](http://cvxr.com/cvx/)
- Gurobi Solver [http://www.gurobi.com/](http://www.gurobi.com/)

## Installation Instructions

### SReachTools

1. Clone the SReachTools repository (or download the latest zip file from
   [Releases](https://github.com/unm-hscl/SReachTools/releases))

### Model Parametrix Toolbox (MPT)

1. Install MPT according to the [installation instructions](https://www.mpt3.org/Main/Installation)

### CVX

1. Download CVX "Standard bundles, including Gurobi and/or MOSEK" for the 
   appropriate operating system. [Download page](http://cvxr.com/cvx/download/)
1. Follow the remaining [installation instructions](http://cvxr.com/cvx/doc/install.html)

### Gurobi

1. [Dowload Gurobi 7.5.2](http://www.gurobi.com/downloads/gurobi-optimizer)
1. Obtain a [free academic license](https://user.gurobi.com/download/licenses/free-academic)
1. Follow the remaing [instructions](http://www.gurobi.com/documentation/8.0/quickstart_mac/matlab_setting_up_gurobi_f.html) for setting up gurobi with MATLAB

## Instructions for Repeating Results

### Initialization

After launching an instance of MATLAB follow these steps:

1. Add MPT to the path by running `tbxmanager restorepath;`
1. Navigate to the CVX directory and run the following commands to setup CVX
   ```
   cvx_solver Gurobi;
   cvx_saveprefs;
   cvx_setup;
   ```
1. Navigate to the SReachTools directory and run `srtinit;` to initialize it

### Reproducing Figures

- Navigate to the directory of the repeatability codes
- Run FigureX.m to generate Figure X in the paper. X = {2,3,4,5}

## Contact details

* Abraham P. Vinod ([aby.vinod@gmail.com](mailto:aby.vinod@gmail.com))
* Joseph D. Gleason ([gleasonj@unm.edu](mailto:gleasonj@unm.edu))
* Meeko M. K. Oishi ([oishi@unm.edu](mailto:oishi@unm.edu))

