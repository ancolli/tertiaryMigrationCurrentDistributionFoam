# v1.0.0

Cite as: Ind. Eng. Chem. Res. 60 (2021) 11927
DOI: https://doi.org/10.1021/acs.iecr.1c01884

# tertiaryMigrationCurrentDistributionFoam
The process of simulating tertiary current distribution, including the effect of ion migration, in electrochemical reactors is described using OpenFOAM, the custom-developed solver potentialMigrConcentrationFoam, and the codedMixed boundary condition. The workflow demonstrates how to pre-process, run, and post-process a basic case in a 2D domain.

The proposed strategy supports the input of kk-th electrochemical reactions per electrode and is fully functional in OpenFOAM version 8 and version 2312.

# Disclaimer
This offering is not approved or endorsed by OpeFOAM Foundation, producer and distributor of the OpenFOAM software via www.openfoam.org.

# Usage
In applications (A) you will find the scripts to compile the solver and post-processing utilities in order to solve tertiary current distribution including the migration effect.
In tutorial (B) you will find an example of a 2D parallel-plate cell. 

Valid for k-th species without the need of modifying the solver (Authomatic generation of species from speciesProperties dictionary).  

# #  A) Applications
**1.**  Solver  
_A)_ Paste applications/utilities/solvers/potentialMigrConcentrationFoam inside OpenFOAM user directory (Applications/Utilities/Solvers).  
_B)_ Open a terminal inside potentialMigrConcentrationFoam.  
_C)_ Run wmake.  
**2.**  Current density due to a given species  
_A)_ Paste applications/utilities/postProcessing/speciesMigrationWallFlux inside OpenFOAM user directory (Applications/Utilities/postProcessing).  
_B)_ Open a terminal inside speciesMigrationWallFlux.  
_C)_ Run wmake.  
**3.**  Total current density  
_A)_ Paste applications/utilities/postProcessing/migrationWallFlux inside OpenFOAM user directory (Applications/Utilities/postProcessing).  
_B)_ Open a terminal inside migrationWallFlux.  
_C)_ Run wmake.  


# #  B) Tutorial
**1-** Paste tutorial inside OpenFOAM user directory (Run/Tutorials).  
**2-** Enter to tutorial and open a Terminal.  
**3-** Modify transport properties (conductivity, kinetic parameters, concentrations, equilibrium potentials, Sct, diffusion coefficients) inside constant/transportProperties.  
**4-** Modify control properties (kind of electric control, under-over relaxation parameter, maxPot, tolerances) inside constant/controlProperties.   
**5-** Run ./Allrun.    

