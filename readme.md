# 2d tissue simulator

author: Attila Gabor  
contact: attila.gabor@uni-heidelberg.de  
date: 2022-04-04


This repo contains the code and data used to generate the insilico data for the 
paper:  
Tanevski et al (2022) **Explainable multi-view framework for dissecting inter-cellular signaling from highly multiplexed spatial data**, *MSB*, in print. 


# Dependencies
Note that this repo is not an R package. Dependencies are scattered in scripts 
and *not* installed automatically, but the basic dependencies are the following

```{r}
library(tidyverse)   # general data manipulation and plotting
library(igraph)      # used for generating the intracellular networks of cell types
library(CellNOptR)   # (from Bioconductor) used for generating the intracellular networks of cell types
library(R6)          # for the object oriented core
```

# 1. Basic idea

The simulator solves a reaction diffusion system on a 2D lattice.
Each grid point can be occupied by a cell type or left empty. If a node is occupied
by a cell type, then the cell-type specific rules describes how the species
are produced, how they interact with each other and how they degrade. On the empty grid points
the species are only diffusing and degrading. 

# 2. The simulator

At the heart of the simulator is an IMC object ('./IMC/insilico_data_generation/IMC_model.R'). 
The IMC object contains a 2D lattice filled with cells (nodes). A node can be empty or 
filled with a cell type. This object is responsable for initialization, running the simulation,
plotting the results (time course and 2D density plots) and exporting the states. 

The IMC object is initialized with a segmentation file and a config file.  
The IMC object can then solve the reaction-diffusion system
parameterized in the config file. It can plot the distribution of compounds or
export the data at certain time points. 

## 2.1 Segmentation file

A segmentation file contains the layout of the cell types in a plain csv file.  
Example: 'IMC/segmentation/Tissues/Tissue1_grid_preproc.csv'.  It contains cell types and 
empty cells. Note that the user has to configure the reaction kinetics for each
node types: all cell types + empty node.

## 2.2 The config file

The config file describes how the simulator should work.
An example can be found at "/IMC/insilico_data_generation/general_config_v8_correctIC.R". 
The config file has to declare the following variables:  

- `diffCoeff`: vector of diffusion coefficients for each modeled species
- `initState`: list (length = n_node_type) initial state for each species in each node type
- `degradCoeff`: list (length = n_node_type) of vectors storing the degradation coefficient of each 
species
- `prodCoeff`: list (length = n_node_type) of matrices. Describe how species influence the activation/production of each other. 


### 2.2.1 So what is the first 242 lines for? 

For the simulator, we had to define cell types with different intracellular signaling pathways,
different set of receptors and ligand-production profiles. 

To do this, we build a general protein-protein interaction network using a random 
graph. Then subsets of this graph serve as cell-type specific models. 
Finally, using some heuristics (e.g. the nodes which have no inputs are ligands),
we define ligands, receptors and transcription factors, which controls the ligand 
production. 



# 4. Main files and scripts: 

- `case_studies/insilico_IMC_generator_TIGER_config_v8_single.R` driver script to run the simulator. Runs one layout.
- `case_studies/insilico_IMC_generator_TIGER_config_v8_batch.R` driver script to run the simulator. Runs all 4 layouts in parallel. 
- `./IMC/insilico_data_generation/IMC_model.R` : the core of the simulator
- `./IMC/segmentation/Tissues/Tissue1_grid_preproc.csv` example segmentation file, contains
the layout of cell types and inter cellular space. 
- `/IMC/insilico_data_generation/general_config_v8_correctIC.R` config file for 
the definition of the reaction-diffusion system. 


# 5. How to start

1. clone/download the repo
2. install the dependencies
3. run the case study script line by line
`case_studies/insilico_IMC_generator_TIGER_config_v8_single.R`

license: GPL-v3

