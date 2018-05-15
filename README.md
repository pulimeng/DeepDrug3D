# DeepDrug3D

DeepDrug3D is a tool to classify the protein pocket to be ATP-binding or Heme-binding given the binding residue numbers and the protein structure.

This README file is written by Limeng Pu

![eg_image](https://github.com/pulimeng/DeepDrug3D/blob/master/image/eg_grid.png)

An example of binding grid generated, pdb ID: 1a2sA.

# Prerequisites

1. Python 2.7+
2. numpy 7.8.2 or higher
3. scipy 0.13.3 or higher
4. scikit-learn 0.19.0 or higher
5. Openbabel 2.3.1 or higher (if you are using Anaconda, can be installed with `conda install -c openbabel openbabel`)
6. tensorflow-gpu 
8. CUDA 7.5 or higher
9. keras 2.1.4 or higher

For the installation instruction please refer to the corresponding project site.

# Usage

The package provides both prediction and training modules. 

1. The prediction module 
It uses the pdb file and an auxilary input file, which contains biniding residue numbers and center of the ligand/pocket, as input files. The center in the auxilary input file is not necessary. If the center is not provided, the model will calculate the pocket center and use it as the ligand center. An example of the auxilary file is provided in `example_aux.txt`. The prediction modle is available at `TODO insert model website here!`
To use the prediction module, run `python predict.py --protein your_protein.pdb --aux your_auxilary_file.txt`.
  - The `--protein` contains the full path to the pdb file you wish to classify.
  - The `--aux` is the auxilary file with binding residue numbers and center of ligand (optional).
  - Two files will be generated along the process, namely `your_protein_trans.pdb` and `your_protein_trans.mol2` under the current working directory. These files are the transformed (moved to the provided center and aligned with the principal axes of the pocket) protein. They will be used during the preceding processes. If you do not wish to keep them, you can just delete them after the getting the results.
  - The output will be printed as two probabilities that how likely the given pocket is to bind with ATP and Heme.
  - The entire process may take upto 30 minutes to finish since the grid point generation is very time consuming.
  - The DFIRE potentials calculation uses the module provided by `A Knowledge-Based Energy Function for Protein−Ligand, Protein−Protein, and Protein−DNA Complexes by Zhang et al.` since it is written in Fortran, which is faster than our own implementation in Python.
  
2. The training module
The 
