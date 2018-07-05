def ramachandran(u,start,finish,r):
    import MDAnalysis as mda
    import numpy as np

    # Creates an array of resids that are used as indicies
    resids = r.residues.resids

    # Creates an array that gives atom groups for respective angles
    dh_atoms = []
    for resid in resids:
        resname = u.residues[resid-1]
        dh_atoms.append(resname)
    dh_atoms = np.array(dh_atoms)

    # Creates an array of time steps that contain [phi,psi] for each residue
    dh_angles = []
    for ts in u.trajectory[start:finish]:
        step = []
        for resid in resids:
            if resid > 1 and resid < u.atoms.select_atoms("protein").n_residues:
                bb = u.atoms.select_atoms("resid {}:{} and name N CA C".format((resid-1),(resid+1)))
                phi = bb.atoms[2:6].dihedral.value()
                psi = bb.atoms[3:7].dihedral.value()
                step.append((phi, psi))
            else:
                pass
        dh_angles.append(step)
    dh_angles = np.array(dh_angles)

    # Reference atom groups with phi_psi(u)[0][residue] --> [<Atom 23: C ...> <Atom 24: N ...> ...]
    # Reference angles with phi_psi(u)[1][time_step][residue] --> [phi,psi]
    return dh_atoms, dh_angles
