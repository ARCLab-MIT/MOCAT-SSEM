# Stable and Sustainable Orbital Capacity Solutions in Low Earth Orbit

Note: the files will be uploaded after the review and approval from the Department of the Air Force Artificial Intelligence Accelerator. 

## Description

The main contributions of this work are the following: 
(1) An improved MOCAT-SSEM that incorporates a novel debris spreading feature is presented. 
(2) Debris tend to accumulate in high-altitude orbits. 
(3) A novel optimization capability is introduced in MOCAT-SSEM to compute the LEO optimal orbital carrying capacity solutions, considering stable equilibrium points and sustainability constraints. 
(4) A first-ever comparison and validation analysis is performed between a SSEM and MC model, specifically MOCAT-SSEM and MOCAT-MC. 
(5) The equilibrium solutions, which are retrieved in MOCAT-SSEM, exist in MOCAT-MC. 
(6) The optimal orbital carrying capacity solutions, which are retrieved in MOCAT-SSEM, can be reproduced in MOCAT-MC.
A three-species MOCAT-SSEM is considered, with active satellites (S), derelicts (D), and debris (N).


## System requirements

Matlab with Optimization Toolbox.


## Installation guide

Download this folder and run main file: "optimal_capacity.m" to compute the LEO optimal orbital carrying capacity solutions, considering stable equilibrium points and sustainability constraints.


## Demo/Output with a three-species MOCAT-SSEM

The main file "optimal_capacity.m" is set to reproduce the results of Case 2 and 3 (by commenting in/out the lines of code referred to the Inequality constraints in PrcAll.m) reported in Table 2 and Figure 2. 
Run time is usually in the order of seconds. 


## Validation with MOCAT-MC

The folder "setup_files_MC" contains the MOCAT-MC configuration files to recreate the simulation scenarios, by using MOCAT-MC which is openly accessible [here](https://github.com/ARCLab-MIT/MOCAT-MC).
The folder "1_to_run_MC" contains the scripts to run MOCAT-MC. The folder "2_to_create_summary" contains the scripts to create a summary file of the MOCAT-MC output for comparison with MOCAT-SSEM.


## Citing

If you find this project research useful, please cite our work:

```
@article{lavezzi_acta_25,
author = {Lavezzi, Giovanni and Jang, Dan and Linares, Richard},
journal={Under peer review at Acta Astronautica}, 
title = {Stable and Sustainable Orbital Capacity Solutions in Low Earth Orbit},
year = {2025},
volume = {},
number = {},
pages = {},
}
```

<!-- doi = {10.21203/rs.3.rs-4927491/v1}, -->
<!-- Check our paper, currently under peer review in Acta Astronautica. Link to the [preprint](https://www.researchgate.net/publication/383876801_A_Stable_and_Sustainable_Maximum_Orbital_Capacity_in_Low_Earth_Orbit). -->


## Acknowledgments

Research was sponsored by the Department of the Air Force Artificial Intelligence Accelerator and was accomplished under Cooperative Agreement Number FA8750-19-2-1000. The views and conclusions contained in this document are those of the authors and should not be interpreted as representing the official policies, either expressed or implied, of the Department of the Air Force or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for Government purposes notwithstanding any copyright notation herein.
The authors acknowledge the MIT SuperCloud and Lincoln Laboratory Supercomputing Center for providing HPC resources that have contributed to the research results reported within this paper.

