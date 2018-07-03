# Calculates dihedral angles of residues specified by user
# Requires a usable atom group
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np

# Test topology and trajectory
u = mda.Universe('md_0_1.gro', 'md_0_1_noPBC1.xtc')
bb = u.select_atoms("name CA C N")

# Creates an array of atom groups that correspond to the angles
dh_atoms = []
for i in range(0,len(bb.atoms)-4,3):
    dh_atoms.append(bb.atoms[i:(i+5)])
dh_atoms = np.array(dh_atoms)

# Creates a numpy array of dihedral angles
# [ [[phi,psi]  [[phi,psi]   ...   [[phi,psi]        resid = 1
#    [phi,psi]   [phi,psi]   ...    [phi,psi]        resid = 2
#       ...         ...      ...       ...              ...
#    [phi,psi]]  [phi,psi]]  ...    [phi,psi]] ]     resid = n_residues
#     ts = 0      ts = 1     ...   ts = n_frames
dha = []
for ts in u.trajectory:
    step = []
    for i in range(2,len(bb.atoms)-4,3):
        phi = bb.atoms[i:(i+4)].dihedral.value()
        psi = bb.atoms[(i+1):(i+5)].dihedral.value()
        step.append((phi, psi))
    dha.append(step)
dha = np.array(dha)

fig = plt.figure(figsize=(10,10))

for i in range(0,len(u.trajectory)):
    plt.plot(dha[i, :,0], dha[i, :,1], 'bx', label="Ramachandran Plot")

plt.axis([-180,180,-180,180])
plt.axhline(0, color='k', lw=1)
plt.axvline(0, color='k', lw=1)
plt.xticks(np.arange(-180,181,60))
plt.yticks(np.arange(-180,181,60))
plt.title("Ramachandran Plot of Lysozyme (PBD: 1AKI)")
plt.xlabel(r"$\phi$ (deg)")
plt.ylabel(r"$\psi$ (deg)")
plt.show()
