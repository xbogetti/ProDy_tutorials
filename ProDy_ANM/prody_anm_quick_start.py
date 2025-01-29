from prody import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size':10})

#Load structure. In this example, we fetch an example PDB file: 4ake.pdb. You may want to change the name of the file that you want to analyze.
open_aa = parsePDB('4ake', compressed=False)
#Select the alpha-carbons of chain A
open_ca = open_aa.select('calpha and chain A')

#Initiate a ANM object, build the connectivity matrice using the buildHessian() function, and calculate the ANM modes. Here, we calculate the 20 slowest ANM modes (the default number of modes to be calculated). You may change the number of modes to be calculated. Refer to the notebook tutorials for more details.
anm_open = ANM('open_AKE')
anm_open.buildHessian(open_ca)
anm_open.calcModes()

# The following commands can be used to extend the anm modes from alpha-carbons to all atoms. This step may be needed for visualization using softwares such as PyMOL and VMD.
aa_anm, aa_atoms = extendModel(anm_open, open_ca, open_aa)
## The follwing command cn be used to save an NMD file, which is a plain text files that contain the normal modes and system coordinate data. NMD files can be visualized using Normal Mode Wizard within the VMD software. ProDy functions writeNMD() and parseNMD() can be used to read and write NMD files.
writeNMD('4ake_CA_anm', anm_open, open_ca)
# Or you can save the NMD file for the extended ANM modes:
writeNMD('4ake_all_atom_anm', aa_anm, aa_atoms)

##Save the ANM modes for plotting later. If you would like to save the ANM modes calculated above, please use the following commands:
saveAtoms(open_ca, '4ake_CA_atoms')
saveModel(anm_open, '4ake_CA')
##Save the extended ANM modes:
saveAtoms(aa_atoms, '4ake_all_atoms')
saveModel(aa_anm, '4ake_all')

##If you have already calculated the ANM modes of a molecule and saved the model, you can load the modes and plot them:
mol_ca = loadAtoms('4ake_CA_atoms')
anm_open = loadModel('4ake_CA.anm.npz')

##Plotting the Mean Squared Fluctuations (MSF) using your personalized plotting scripts. You may use the following codes to plot the MSF of specific ANM mode(s). In this following example, we plot the MSF of ANM modes 1 to 5 and the MSF can be obtained by the function: .
msfs = calcSqFlucts(anm_open)
resnum = np.arange(0, msfs.shape[0], 1)
##This PDB starts with residue 1
resnum = resnum + 1

plt.figure(figsize=(10, 4))
plt.plot(resnum, msfs, c = 'k', linewidth = 2)
#plt.xlim(0, 355)
plt.xlabel('Residue number')
plt.ylabel('MSF')
plt.tight_layout()
plt.savefig('4ake_msfs.png', dpi=300)

