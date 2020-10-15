# DeepDrug3D

DeepDrug3D is a tool to predict the protein pocket to be ATP/Heme/other-binding given the binding residue numbers and the protein structure.

If you find this tool useful, please star this repo and cite our paper :)

Pu L, Govindaraj RG, Lemoine JM, Wu HC, Brylinski M (2019) DeepDrug3D: Classification of ligand-binding pockets in proteins with a convolutional neural network. PLOS Computational Biology 15(2): e1006718. https://doi.org/10.1371/journal.pcbi.1006718

This README file is written by Limeng Pu. **This is a newer version of the implmentation. A intermediate results visualization function has been added. And the model has been retrained using Pytorch instead of Keras.**

<p align="center">
    <img width="400" height="400" src="./image/1a2sA.png">
</p>

An example of binding grid generated, pdb ID: 1a2sA, atom type: C.ar. Red --> low potentials while Blue --> high potentials.

# Prerequisites
1. System requirement: Linux (DFIRE potential calculation only runs on Linux. Tested on <em>Red Hat Enterprise Linux 6</em>)
2. The dependencies are provided in `environment.yml`. Please change line 9 and line 14 in the file according to your system. To install all the dependencies run `conda env create -f environment.yml`.
3. In case of failuer to run DFIRE or Pytorch, separate the data-generation process and trainig/prediction process. Make two distinct environment for each of them. Then store the generated data for training/prediction.

# Usage

The package provides data-generation, prediction, and training modules.

1. Data generation

This step generates data for training/prediction while provides intermediate results for visualization.

The input contains a .pdb file of the protein and auxilary file marking the binding residue numbers and center of the pocket (optional). An example of the auxilary file is provided in `example_aux.txt`.

2. Prediction

It uses the pdb file and an auxilary input file, which contains biniding residue numbers and center of the ligand/pocket, as input files. The center in the auxilary input file is not necessary. If the center is not provided, the model will calculate the pocket center and use it as the ligand center.  The trained model is available at `https://osf.io/enz69/`
To use the prediction module, run `python predict.py --protein your_protein.pdb --aux your_auxilary_file.txt --r 15 --N 31`.
  - `--protein` contains the full path to the pdb file you wish to classify.
  - `--aux` is the auxilary file with binding residue numbers and center of ligand (optional).
  - `--r` and `--N` are the radius of the grid and number of points along the dimension of the grid. The default settings are r = 15 and N = 31. This setting yeilds a 32 x 32 x 32 grid. This can be changed by setting r and N.
  - Two files will be generated along the process, namely `your_protein_trans.pdb` and `your_protein_trans.mol2` under the current working directory. These files are the transformed (moved to the provided center and aligned with the principal axes of the pocket) protein. They will be used during the later processes. If you do not wish to keep them, you can just delete them after the getting the results.
  - The output will be printed as three probabilities that each represents the likelihood of the pocket being an ATP/Heme/other binding pocket.
  - The entire process may take upto 30 minutes to finish since the grid point generation (mostly the potential calculation) is very time consuming.
  - The DFIRE potentials calculation uses the module provided by `A Knowledge-Based Energy Function for Protein−Ligand, Protein−Protein, and Protein−DNA Complexes by Zhang et al.` since it is written in Fortran, which is faster than our own implementation in Python.
  
3. Training

In order to use our model to train your own dataset, you have to conert your dataset, which will be pdbs to voxel representation of protein-ligand biniding site. The trainig module can be runned as `python train.py --alist deepdrug3d_atp.lst --hlist deepdrug3d_heme.lst --vfolder deepdrug3d_voxel_data --bs batch_size --lr inital_learning_rate --epoch number_of_epoches --output deepdrug3d`.
  - `--alist` is the list of the full path to the ATP binding voxel data while `--hlist` is the list of the full path to the Heme binidng voxel data.
  - `--vfolder` is the folder contains all the voxel data, which contains numpy array (.npy) for each protein-ligand pair.
  - `--bs`, `--lr`, `--epoch` is the hyperparameters related to the model. Recommanded values are 64, 0.00001, 30.
  - If no output location is provided, the model will be saved to the current workding direcotry as 'deepdrug3d.h5'.
  
# Dataset

We provided our dataset we used for the training at `https://osf.io/enz69/`, which are the voxel representations of ATP, Heme, and other along with the class label file.
