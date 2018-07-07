def dihedral_calc(ags, start=None, stop=None, step=None):
    """Calculates phi and psi angles for a list of AtomGroups over trajectory.

    Parameters
    ----------
    ags : AtomGroups
        must be a list of one or more AtomGroups containing 5 atoms in the
        correct order (i.e. C-N-CA-C-N)
     start : int, optional
         starting frame of analysis
     stop : int, optional
         last frame of analysis plus 1
     step : int, optional
         step size between frames for analysis

    Returns
    -------
    angles : numpy.ndarray
        An array of time steps which contain (phi,psi) for all AtomGroups.
    """

    phi_sel = [ag[:4].dihedral for ag in ags]
    psi_sel = [ag[1:].dihedral for ag in ags]

    angles = np.array([[(phi.value(), psi.value())
             for phi,psi in zip(phi_sel,psi_sel)]
             for ts in ags.universe.trajectory[start:stop:step]])

    return angles


def ramachandran(r, start=None, stop=None, step=None):
    """Generates time series of phi and psi angles for all residues specified.

    Parameters
    ----------
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

    Time steps can be found with x[1][ts] and angles for a specific residue
    can be found with x[1][ts][resid+1]
    """

    protein = r.universe.atoms.select_atoms("protein")
    resids = r.residues.resids
    bb_sel = [protein.select_atoms(
                  "name C and resid {}".format(resid-1),
                  "name N and resid {}".format(resid),
                  "name CA and resid {}".format(resid),
                  "name C and resid {}".format(resid),
                  "name N and resid {}".format(resid+1)
                  )
              for resid in resids
              if 1 < resid < len(protein.residues)]
    angles = dihedral_calc(bb_sel,start=start,stop=stop,step=step)

    return r.residues, angles
