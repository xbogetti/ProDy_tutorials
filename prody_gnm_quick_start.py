from prody import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size':10})

#Load structure. In this example, we load the chain A.
pdb = parsePDB('1p38.pdb', chain='A')
mol = pdb.protein 
#Select the alpha-carbons
mol_ca = mol.ca

#Initiate a GNM object, build the connectivity matrice using the buildKirchhoff() function, and calculate the GNM modes. Here, we calculate the 10 slowest GNM modes. You may change the number of modes to be calculated. Refer to the notebook tutorials for more details.
gnm = GNM('example')
gnm.buildKirchhoff(mol_ca)
gnm.calcModes(n_modes=10)

# The following commands can be used to extend the gnm modes from alpha-carbons to all atoms. This step may be needed for visualization using softwares such as PyMOL and VMD.
gnm_aa, atoms_all = extendModel(gnm, mol_ca, mol)
# The follwing command cn be used to save an NMD file, which is a plain text files that contain at least normal mode and system coordinate data. NMD files can be visualized using Normal Mode Wizard within the VMD software. ProDy functions writeNMD() and parseNMD() can be used to read and write NMD files.
writeNMD('1p38_gnm_aa.nmd', gnm_aa, mol)

#Save the GNM modes for plotting later. If you would like to save the GNM modes calculated above, please use the following commands:
saveAtoms(mol_ca, '1p38_CA_atoms')
saveModel(gnm, '1p38_CA')
#Save the GNM modes for the extended GNM modes:
saveAtoms(mol, '1p38_all_atoms')
saveModel(gnm_aa, '1p38_all')

##If you have already calculated the GNM modes of a molecule and saved the model, you can load the modes and plot them:
mol_ca = loadAtoms('1p38_CA_atoms')
gnm = loadModel('1p38_CA.gnm.npz')

##Plotting each mode using your personalized plotting scripts. You may use the following codes to plot the mode shapes (eigenvectors of a specific mode). In this following example, we plot eigenvectors of GNM mode 1 and the eigenvectors can be obtained by the function: getEigvecs().
mode = gnm[0].getEigvecs()
resnum = np.arange(0, mode.shape[0], 1)
##This PDB starts with residue 4
resnum = resnum + 4

plt.figure(figsize=(10, 4))
plt.plot(resnum, mode, c = 'k', linewidth = 2)
plt.axhline(y=0, c = 'darkgrey', linestyle = '--', linewidth = 2)
plt.xlim(0, 355)
plt.xlabel('Residue number')
plt.ylabel('GNM mode 1 profile')
plt.tight_layout()
plt.savefig('1p38_gnm_mode1.png', dpi=300)

