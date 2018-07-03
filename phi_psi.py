def phi_psi(u):
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np

    # Selects relevent backbone atoms
    bb = u.select_atoms("name CA C N and not ((name N CA and resid 1) or (name CA C and resid {}))".format(len(u.residues)))

    # Creates an array that gives atom groups for respective angles
    dh_atoms = []
    for i in range(0,len(bb.atoms)-4,3):
        dh_atoms.append(bb.atoms[i:(i+5)])
    dh_atoms = np.array(dh_atoms)

    # Creates an array of time steps that contain [phi,psi] for each residue
    dh_angles = []
    for ts in u.trajectory:
        step = []
        for i in range(0,len(bb.atoms)-4,3):
            phi = bb.atoms[i:(i+4)].dihedral.value()
            psi = bb.atoms[(i+1):(i+5)].dihedral.value()
            step.append((phi, psi))
        dh_angles.append(step)
    dh_angles = np.array(dh_angles)

    # Reference atom groups with phi_psi(u)[0][residue] --> [<Atom 23: C ...> <Atom 24: N ...> ...]
    # Reference angles with phi_psi(u)[1][time_step][residue] --> [phi,psi]
    return dh_atoms, dh_angles
