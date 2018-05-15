# DeepDrug3D

DeepDrug3D is a tool to classify the protein pocket to be ATP-binding or Heme-binding given the binding residue numbers and the protein structure.

This README file is written by Limeng Pu

# Prerequisites

1. Python 2.7+
2. numpy 7.8.2 or higher
3. scipy 0.13.3 or higher
4. scikit-learn 0.19.0 or higher
5. Openbabel 2.3.1 (if you are using Anaconda, can be installed with `conda install -c openbabel openbabel`)
6. tensorflow-gpu 
8. CUDA 7.5 or higher
9. keras 2.1.4 or higher

# Usage

The package provides both prediction and training modules. The prediction module uses the pdb file and an auxilary input file, which contains biniding residue numbers and center of the ligand/pocket. The center is not necessary. If the center is not provided, the model will calculate the pocket center and use it as the ligand center.
