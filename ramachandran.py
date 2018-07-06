def ramachandran(u,r, start=None, stop=None, step=None):
    """Generates a time series of phi and psi angles for all residues specified.

    Parameters
    ----------
    u : Universe object
        because of how dihedrals are calculated, residues that may not be
        specified must be able to be referred to
    r : AtomGroup or ResidueGroup or SegmentGroup
        some group of interest that a list of residues can be generated from
    start : int, optional
        starting frame of analysis
    stop : int, optional
        last frame of analysis plus 1
    step : int, optional
        step size between frames for analysis

    Returns
    -------
    r.residues : ResidueGroup
        List of residues that correspond to the list of angles
    angles : numpy.ndarray
        An array of time steps which contain (phi,psi) for all residues in r.

    Notes
    -----
    The list of residues can be found with x[0] (i.e. <ResidueGroup with n
    residues>) and specific residues can be foudn with x[0][resid+1].

    Time steps can be found with x[1][ts] and angles for a specific residue can
    be found with x[1][ts][resid+1]
    """

    import MDAnalysis as mda
    import numpy as np

    resids = r.residues.resids
    bb_sel = [u.atoms.select_atoms(
                  "name C and resid {}".format(resid-1),
                  "name N and resid {}".format(resid),
                  "name CA and resid {}".format(resid),
                  "name C and resid {}".format(resid),
                  "name N and resid {}".format(resid+1)
                  )
              for resid in resids
              if (resid > 1 and
                  resid < len(u.atoms.select_atoms("protein").residues))]

    phi_sel = [bb_sel[i][:4].dihedral for i in range(len(bb_sel))]
    psi_sel = [bb_sel[i][1:].dihedral for i in range(len(bb_sel))]

    angles = np.array([[(phi_sel[i].value(), psi_sel[i].value())
             for i in range(len(phi_sel))]
             for ts in u.trajectory[start:stop:step]])

    return r.residues, angles
