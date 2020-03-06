# Thesis_MarcLatour
This repository contains all the relevant files used throughout this thesis.

## Prerequisites for running everything on this repo
- Perform an AD build of SU2_feature_turbo_bodyforce from ``https://github.com/su2code/SU2/tree/feature_turbo_bodyforce`` . This includes Hall's body force in the direct and adjoint solvers and will allow all simulations to be run.
- Download and compile UMG2 from the turbomachinery BitBucket. Ask Matteo Pini for access if you do not have so yet.

## Overview of contents of this repo

#### Meshing

This folder contains all of the meshes that were used. It includes structured and unstructured versions of the body force simulation mesh as well as the unstructured mesh used for the Whittle simulation. Mesh creation folders are included as well, meaing that once UMG2 is installed these meshes can be modified and created again.

#### Direct Implementation

Everything related to the BFM in the direct solver is included in this folder. Camber line generation and airfoil coordinates are included as well, which allows both to be made for use in body force and physical blade simulations respectively. The remaining folders contain all the test cases that are presented in the report. They include all the simulations that were run, including configuration files, meshes, and all output.

#### Adjoint Implementation

Three test cases are presented in this folder, one for each angle of attack. In each of these folders the results of simulations for non-registered body forces and registered body forces are presented. For each of these the objective function is given in the square brackets in the folder name. When DRAG is used as objective function, then it uses the self-defined cumulative velocity in the y-direction (see CEulerSolver::PressureForces in the solver_direct_mean.cpp file).

#### Verification

All files used for analysis of teh results for both the direct and adjoint simulations is in this folder. For each simulation from the direct and adjoint solvers the output is included in these folders as .dat and .png files. For each test case a python file is included that contains the code used to create the analyses and images used in the report. Meridional flow files and history files were taken from the respective test cases.

#### Resources

Included in this folder is the implementation of Hall's body force model in openFOAM (obtained from Jeff Defoe's team at the University of Windsor) as well as a couple of macros used in Tecplot to create output for physical blade simulations. It also includes the virids colormap that can be used in Tecplot, this is a colormap that can be discerned in black and white as well. To support research I have included all the papers, literature, and other sources used during the literature study and thesis. Finally, there is a folder with all the images present in the thesis.
