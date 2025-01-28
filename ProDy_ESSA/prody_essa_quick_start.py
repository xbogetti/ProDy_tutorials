from prody import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt ##These three packages are for plotting

matplotlib.rcParams.update({'font.size':10})

#Load structure. In this example, we fetch the structure 1pzo.pdb from the PDB database.
pdb = parsePDB('1pzo', compressed=False)

#Initiate an ESSA object
essa = ESSA()

#Set the system.
#In the structure 1pzo, there are two identical allosteric inhibitors (residues 300 and 301 of chain A). ESSA can be informed about them when setting the system as follows:
essa.setSystem(pdb, lig='A 300 A 301')
#Or if you don't have a ligand, you can use the command: 
#essa.setSystem(pdb)
#If you don't have enough memory, lowmem=True should be chosen, thus only storing eigenvalues/eigenvectors for each perturbed model:
#essa.setSystem(pdb, lowmem=True)

#Scanning can be done by scanResidue() function
essa.scanResidues()

#Save the z-scores of each residue as a text file
essa_zscores = essa.getESSAZscores()
writeArray('essa_zscores_1pzo.txt', essa_zscores, format='%.5f')
#Alternatively, if you would like to save the z-scores as an numpy array (.npy) file, the following command can be used:
essa.saveESSAZscores()

#Save z-scores to the b-factor column of the PDB file for visulization:
essa.writeESSAZscoresToPDB()


#Quick plotting
resnum = np.arange(0, essa_zscores.shape[0], 1)
##This PDB starts with residue 26
resnum = resnum + 26

plt.figure(figsize=(10, 4))
plt.plot(resnum, essa_zscores, c = 'k', linewidth = 2)
plt.xlim(26, 300)
plt.xlabel('Residue number')
plt.ylabel('ESSA z-scores')
plt.tight_layout()
plt.savefig('1pzo_essa_zscores.png', dpi=300)
