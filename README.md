# MOCAT-SSEM: MIT Orbital Capacity Assessment Tool - Source-Sink Evolutionary Model

# Beta - Version 0.1

MOCAT-SSEM investigates the evolution of the space objects population in LEO by exploiting a new probabilistic source-sink model, with the objective of estimating the LEO orbital capacity. This is carried out through the long-term propagation of the proposed source-sink model, which globally takes into account different object species, such as active satellites, derelict satellites, debris, and additional subgroups. Since the SOs are propagated as species, the information about single objects is missing, but it allows to be computationally fast and provide essential information about projected future distribution of SOs in the space environment for long prediction horizon. MOCAT-SSEM is developed in MATLAB.

<!--
This repository provides a comprehensive package for the replication of different versions of MOCAT-SSEM using MATLAB, described in the following papers: (TO BE ADDED).
-->

## Structure repository

### Main Directory

The folder Demo Workbooks contains the following main scripts:

- **`MOCAT_3_Workbook.mlx`**: This Matlab live script serves as the main script for running MOCAT-3. The MOCAT-3 model contains the following species: unslotted satellites (Su), derelict (D), and debris (N). 

- **`MOCAT_4S_Workbook.mlx`**: This Matlab live script serves as the main script for running MOCAT-4S. The MOCAT-4S model contains the following species: slotted satellites (S), unslotted satellites (Su), derelict (D), and debris (N).

- **`MOCAT_Multi_Mass_Species_Workbook.mlx`**: This Matlab live script serves as the main script for a general implementation of the MOCAT-SSEM framework, in which multi-mass species can be created.

Each Matlab live script reports: 
- An introduction of the model, with equations and parameters definition. 
- A walkthrough on how to setup and run the script.
- Several plots to highlight the outputs of the model.

### Work in progress

- Optimizer runs demo
- Implementation of endogenous initial population 
- Mass-binned debris  

### Supporting Directories

Each sub-directory contains the sub-modules used by the workbooks. The sub-modules are scripts/model/data related to differenct aspects, e.g., the atmospheric drag, the initial space objects population, the launch traffic, etc.

## Usage

1. Clone the repository.
2. Navigate to the folder Demo Workbooks in the root directory.
3. Open one of the Matlab livescript workbooks listed above.

## Tutorial

The folder Tutorials contains previews of output solutions retrieved by MOCAT-3 and MOCAT-4S using constant launch rate and static exponential or JB2008-based atmospheric density models, and by MOCAT-Multi-mass-species. These can be accessed online without the need to install/download any software.

<!--
The following website contains an example of the MOCAT-SSEM Matlab livescript functionalities: [https://gdmend.github.io/MOCAT/](https://gdmend.github.io/MOCAT/) . (TO BE CHANGED)
-->

## Acknowledgments

MOCAT’s initial development has been supported by the Defense Advanced Research Projects Agency (DARPA) under Grant N66001-20-1-4028, NASA’s Office of Technology and Strategy, under Grant Number NNH22ZDA001N-ESPOD.
