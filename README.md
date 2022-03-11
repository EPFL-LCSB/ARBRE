# ARBRE

The data and scripts contained in this repository allow the user to reproduce the figures 
in the ARBRE manuscript from the original data files
and generate novel pathways predicitons for any of the compounds available in the network.

## Installation

The installation can be completed in less than 10 minutes, including installation of 
dependencies and fetching the data from the git repository.

### Requirements

- python 3
- rdkit environment
- networkx

rdkit is required for balance calculation and visualisation. In case the rdkit is not possible to install,
only the ranked pathways list will be generated.

Please, install rdkit and use rdkit environment for running ARBRE with visualization.
First install anaconda: https://docs.anaconda.com/anaconda/install/index.html

Then install rdkit as described here: https://www.rdkit.org/docs/Install.html.

`$ conda create -c conda-forge -n my-rdkit-env rdkit`
`$ conda activate my-rdkit-env`

Since networkx package is not part of the default rdkit environment, install it to the environment as follows when the environment is activated:

`$ conda install networkx`

### Download repository

`$ git clone https://c4science.ch/source/arbre.git`

If you are installing on macOS, make sure you have Homebrew installed, otherwise you might get "git: 'lfs' is not a git command." error.
Once you installed Homebrew, run

`$ brew install git-lfs`
`$ run git-lfs install`

### Note

Data files are stored using git large file storage (lfs). The make file will install git lfs 
automatically. However, if lfs was not installed previously, the repository has to be 
updated after installation:

`$ git-lfs pull`

This is needed to retrieve the data files from the repository after installation.

# Usage

- create a folder with the name of your project in the "projects" directory (e.g. .../arbre/projects/scoulerine_case1a)
- copy the parameters file from the .../arbre/default_parameters folder to your project folder
- follow the instructions for the parameters.txt adjustment specified in the parameters file (for more details consult the manuscript)

## Pathway search

- Set source and target (first 14 letters of INCHIKEY). Set pathway or network expansion parameter to "p".
- The minimal amount of adjustments for your search is substituting the source and target inchikey by your source and target inchikey.
- In case you are not sure what would be your source compound, create "source_compounds.txt" file that contains 14 letter inchikey for all the potential precursors and set "use_source_compounds_list" parameter to "yes". If the source compounds file is not provided the pathways will be reconstructed towards the Tyr, Phe and Trp.

## Network expansion

- Create "path_compounds.txt" file inside your arbre/projects/{projectname} folder. List the LCSB IDs of the compounds that constitute the pathway one per line (as in morphine example).
- Set pathway or network expansion parameter to "n". Set number of generations around your target pathway.

## For both pathway search and network expansion:

Run the code as following

`$ cd code`

`$ python3 main.py {projectname}`

e.g., 

Searching for the pathways towards several precursors at the same time:

`$ python3 main.py norcoclaurine_precursorsList`

Reproducing case study 1a:

`$ python3 main.py scoulerine_case1a`

Reproducing case study 1b:

`$ python3 main.py benzenecarboperoxic_acid_case1b`

Reproducing case study 2a:

`$ python3 main.py tyrphetrp_case2a`

Reproducing case study 2b:

`$ python3 main.py morphine_case2b`

You will get the following output:

- For pathway search:

ranked pathways list and the visualisation will be in the Output folder

- For network expansion:

1. Output file is located at arbre/output/{projectname}/compound_derivatives.csv, it is ready to be opened with excel. 
2. the .gdf files ready for visualisation in Gephi software are available at arbre/output/{projectname}/visualization/gephifiles

The generation is represented by the color and is labeled in the edges part of the .gdf file in "color VARCHAR" column.

You can install Gephi from https://gephi.org.

# Finishing work with arbre

Deactivate your rdkit environment as follows:

`$ conda deactivate`

# Data updates

The data files used in the repository is the same as the one used in the manuscript, which has been downloaded on 21 June 2021.

For updated network files, please contact the authors of the paper directly.

