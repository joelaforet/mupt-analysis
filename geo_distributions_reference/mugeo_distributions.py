'''MuPT functions for calculating geometric property distributions for 
   CG simulations run in the MuPT ecosystem.'''

__author__ = 'Joe Laforet Jr.'
__email__ = 'jola3134@colorado.edu'

#################################################################################################

from mutils.mutils import get_histogram
import MDAnalysis as mda
import numpy as np
import warnings

from typing import Union
# Bond Distribution
# The following function calculates the distribution of bond lengths in a CG simulation
# for a given set of atoms. User must input which two particle types they are 
# interested in computing the bond data for.

def mupt_bond_distribution(
    gsd_file: str,
    A_name: str,
    B_name: str,
    start: int=0,
    stop: int=-1,
    histogram: bool=False,
    l_min: float=0.0,
    l_max: float=4.0,
    normalize: bool=True,
    bins: Union[float, int, str]="auto",
):
    """Returns the bond length distribution for a given bond pair

    Parameters
    ----------
    gsdfile : str
        Filename of the GSD trajectory.
    A_name, B_name : str
        Name(s) of particles that form the bond pair
        (e.g., "A" and "B" for A-B bonds).
    start : int
        Starting frame index for accumulating bond lengths.
        Negative numbers index from the end. (default 0)
    stop : int
        Final frame index for accumulating bond lengths. (default -1)
    histogram : bool, default=False
        If set to True, places the resulting bonds into a histogram
        and retrums the histogram's bin centers and heights as
        opposed to the actual calcualted bonds.
    l_min : float, default = 0.0
        Sets the minimum bond length to be included in the distribution
    l_max : float, default = 5.0
        Sets the maximum bond length value to be included in the distribution
    normalize : bool, default=False
        If set to True, normalizes the angle distribution by the
        sum of the bin heights, so that the distribution adds up to 1.
    bins : float, int, or str,  default="auto"
        The number of bins to use when finding the distribution
        of bond angles. Using "auto" will set the number of
        bins based on the ideal bin size for the data.
        See the numpy.histogram docs for more details.

    Returns
    -------
    1-D numpy.array  or 2-D numpy.array
        If histogram is False, Array of actual bond angles in degrees
        If histogram is True, returns a 2D array of bin centers and bin heights.

    """
    # Load the GSD file using MDAnalysis
    u = mda.Universe(gsd_file)

    # generate the selection string so we only iterate over the atoms
    # that can possibly participate in the desired bond
    # I.E, if you have particles of type A, B, and C,
    # You don't want to consider bonds that involve B and C atoms

    selection_string = f"name {A_name} or name {B_name}"

    all_bonds = []

    # Select atoms based on the selection string
    for ts in u.trajectory[start:stop]:
        relevent_atoms = u.select_atoms(selection_string)

        # Iterate over all of the bonds in the selection from the universe
        for bond in relevent_atoms.bonds:
            bonded_atoms = list(bond.atoms.names)

            # Check if the bond contains both A_name and B_name
            if A_name in bonded_atoms and B_name in bonded_atoms:
                # if verbose:
                #     #print(f"Bond: {bond}, Atoms: {bond.atoms}, Length: {bond.length()}")
                # else:
                #     #print(f"Bond Length: {bond.length()} between atoms {bond.atoms.ix}")
                all_bonds.append(np.round(bond.length(), 3))

    if histogram:
        if min(all_bonds) < l_min or max(all_bonds) > l_max:
            warnings.warn(
                "There are bond lengths that fall outside of "
                "your set l_min and l_max range. You may want to adjust "
                "this range to include all bond lengths."
            )
        bin_centers, bin_heights = get_histogram(
            data=np.array(all_bonds),
            normalize=normalize,
            bins=bins,
            x_range=(l_min, l_max),
        )
        return np.stack((bin_centers, bin_heights)).T
    else:
        return np.array(all_bonds)

#################################################################################################

# Angle Distribution
# The following function calculates the distribution of bond angles in a CG simulation
# for a given triplet of atoms. User must input which three particle types they are
# interested in computing the angle data for.

from typing import Union

def mupt_angle_distribution(
    gsd_file: str,
    A_name: str,
    B_name: str,
    C_name: str,
    start: int=0,
    stop: int=-1,
    degrees: bool=False,
    histogram: bool=False,
    theta_min: float=0.0,
    theta_max: float=None,
    normalize: bool=False,
    bins: Union[float, int, str]="auto",
):
    """Returns the bond angle distribution for a given triplet of particles

    Parameters
    ----------
    gsdfile : str
        Filename of the GSD trajectory.
    A_name, B_name, C_name : str
        Name(s) of particles that form the angle triplet
        (e.g., "A", "B", and "C" for A-B-C angles).
        They must be given in the same order as they form the angle
    start : int
        Starting frame index for accumulating bond lengths.
        Negative numbers index from the end. (default 0)
    stop : int
        Final frame index for accumulating bond lengths. (default -1)
    degrees : bool, default=False
        If True, the angle values are returned in degrees.
        if False, the angle values are returned in radians.
    histogram : bool, default=False
        If set to True, places the resulting angles into a histogram
        and retrums the histogram's bin centers and heights as
        opposed to the actual calcualted angles.
    theta_min : float, default = 0.0
        Sets the minimum theta value to be included in the distribution
    theta_max : float, default = None
        Sets the maximum theta value to be included in the distribution
        If left as None, then theta_max will be either pi radians or
        180 degrees depending on the value set for the degrees parameter
    normalize : bool, default=False
        If set to True, normalizes the angle distribution by the
        sum of the bin heights, so that the distribution adds up to 1.
    bins : float, int, or str,  default="auto"
        The number of bins to use when finding the distribution
        of bond angles. Using "auto" will set the number of
        bins based on the ideal bin size for the data.
        See the numpy.histogram docs for more details.

    Returns
    -------
    1-D numpy.array  or 2-D numpy.array
        If histogram is False, Array of actual bond angles in degrees
        If histogram is True, returns a 2D array of bin centers and bin heights.

    """

    if not degrees and theta_max is None:
        theta_max = np.pi
    elif degrees and theta_max is None:
        theta_max = 180

    u = mda.Universe(gsd_file)

    # generate the selection string so we only iterate over the atoms
    # that can possibly participate in the desired angle
    # I.E, if you have particles of type A, B, and C,
    # You don't want to consider angles that involve B and C atoms
    # If you are interested in the A-A-A angle

    selection_string = f"name {A_name} or name {B_name} or name {C_name}"

    all_angles = []

    # Select atoms based on the selection string
    for ts in u.trajectory[start:stop]:
        relevent_atoms = u.select_atoms(selection_string)

        # Iterate over all of the angles in the selection from the universe
        for angle in relevent_atoms.angles:
            angled_atoms = list(angle.atoms.names)

            # Only append angles with exact ordering: A_name-B_name-C_name or C_name-B_name-A_name
            if angled_atoms == [A_name, B_name, C_name] or angled_atoms == [C_name, B_name, A_name]:
                # NOTE There appears to be a slight quirk in MDAnalysis where the angle
                # in radians is computed slightly differently than in CMEutils
                # The difference is very small, but it is worth noting
                # See example notebook distributions.ipynb for more details
                if degrees:
                    all_angles.append(np.round(angle.value(), 3))  # by default this is in degrees
                else:
                    all_angles.append(np.round(np.deg2rad(angle.value()), 3))

    if histogram:
            if min(all_angles) < theta_min or max(all_angles) > theta_max:
                warnings.warn(
                    "There are bond angles that fall outside of "
                    "your set theta_min and theta_max range. "
                    "You may want to adjust this range to "
                    "include all bond angles."
                )
            bin_centers, bin_heights = get_histogram(
                data=np.array(all_angles),
                normalize=normalize,
                bins=bins,
                x_range=(theta_min, theta_max),
            )
            return np.stack((bin_centers, bin_heights)).T
    else:
        return np.array(all_angles)
    

#################################################################################################

# Dihedral Distribution
# The following function calculates the distribution of dihedral angles in a CG simulation
# for a given set of four atoms. User must input which four particle types they are
# interested in computing the dihedral data for.

def mupt_dihedral_distribution(
    gsd_file: str,
    A_name: str,
    B_name: str,
    C_name: str,
    D_name: str,
    start: int=0,
    stop: int=-1,
    degrees: bool=False,
    histogram: bool=False,
    normalize: bool=False,
    bins: Union[float, int, str]="auto",
):
    """Returns the dihedral angle distribution for a given set of particles

    Parameters
    ----------
    gsdfile : str
        Filename of the GSD trajectory.
    A_name, B_name, C_name, D_name : str
        Name(s) of particles that form the dihedral angle
        (e.g., "A", "B", "C", and "D" for A-B-C-D dihedrals).
        They must be given in the same order as they form the dihedral
    start : int
        Starting frame index for accumulating bond lengths.
        Negative numbers index from the end. (default 0)
    stop : int
        Final frame index for accumulating bond lengths. (default -1)
    degrees : bool, default=False
        If True, the angle values are returned in degrees.
        if False, the angle values are returned in radians.
    histogram : bool, default=False
        If set to True, places the resulting angles into a histogram
        and retrums the histogram's bin centers and heights as
        opposed to the actual calcualted angles.
    phi_min : float, default = 0.0
        Sets the minimum phi value to be included in the distribution
    phi_max : float, default = None
        Sets the maximum phi value to be included in the distribution
        If left as None, then phi_max will be either pi radians or
        180 degrees depending on the value set for the degrees parameter
    normalize : bool, default=False
        If set to True, normalizes the angle distribution by the
        sum of the bin heights, so that the distribution adds up to 1.
    bins : float, int, or str,  default="auto"
        The number of bins to use when finding the distribution
        of bond angles. Using "auto" will set the number of
        bins based on the ideal bin size for the data.
        See the numpy documentation for more details.
    Returns
    -------
    1-D numpy.array  or 2-D numpy.array
        If histogram is False, Array of actual dihedral angles in degrees
        If histogram is True, returns a 2D array of bin centers and bin heights.

    """

    u = mda.Universe(gsd_file)

    # generate the selection string so we only iterate over the atoms
    # that can possibly participate in the desired dihedral
    # I.E, if you have particles of type A, B, C, and D,
    # You don't want to consider dihedrals that involve B and C atoms
    # If you are interested in the A-A-A-A dihedral

    selection_string = f"name {A_name} or name {B_name} or name {C_name} or name {D_name}"

    all_dihedrals = []

    # Select atoms based on the selection string
    for ts in u.trajectory[start:stop]:
        relevent_atoms = u.select_atoms(selection_string)

        # Iterate over all of the angles in the selection from the universe
        for dihedral in relevent_atoms.dihedrals:
            dihedraled_atoms = list(dihedral.atoms.names)

            # Only append angles with exact ordering: A_name-B_name-C_name or C_name-B_name-A_name
            if dihedraled_atoms == [A_name, B_name, C_name, D_name] or dihedraled_atoms == [D_name, C_name, B_name, A_name]:
                if degrees:
                    all_dihedrals.append(np.round(dihedral.value(), 3))  # by default this is in degrees
                else:
                    all_dihedrals.append(np.round(np.radians(dihedral.value()), 3))

    if histogram:
            bin_centers, bin_heights = get_histogram(
                data=np.array(all_dihedrals),
                normalize=normalize,
                bins=bins,
                x_range=(-np.pi, np.pi),
            )
            return np.stack((bin_centers, bin_heights)).T
    else:
            return np.array(all_dihedrals)
















