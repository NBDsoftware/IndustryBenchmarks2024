import numpy as np
import pathlib
import click
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms


@click.command
@click.option(
    '--prepared',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared PDB file",
)
@click.option(
    '--original',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the original PDB file",
)
def run(prepared, original):
    """
    Compare the prepared PDB with its original form, checking that there
    are not major differences in atom counts and positions.

    Parameters
    ----------
    prepared : pathlib.Path
      A Path to the prepared PDB file.
    original : pathlib.Path
      A Path to the `original` PDB file.
    """
    prep = mda.Universe(prepared)
    prev = mda.Universe(original)
    
    # Align the two structures using just C-alpha atoms
    # align.alignto(prep, prev, select='name CA')

    prep_prot = prep.select_atoms('protein and not resname ACE NME NMA')
    prev_prot = prev.select_atoms('protein and not resname ACE NME NMA')

    for resi, resj in zip(prep_prot.residues, prev_prot.residues):
        print(f"Comparing {resi.resname} {resi.resid} with {resj.resname} {resj.resid}")
        
        if len(resi.atoms) != len(resj.atoms):
            print("WARNING: Number of atoms does not match: ", resi.resname, resi.resid)
            
            # Check if one of the residues has alternative locations
            alt_locs_i = [atom.altLoc for atom in resi.atoms]
            alt_locs_j = [atom.altLoc for atom in resj.atoms]
            
            if any(alt_locs_i):
                print("Residue in prepared structure has alternative locations!")
            if any(alt_locs_j):
                print("Residue in original structure had alternative locations!")
            
        else:
            # Number of atoms match, check if the order of the elements match
            try:
                same_order = all([a1.element == a2.element for a1, a2 in zip(resi.atoms, resj.atoms)])
                if same_order:
                    print("Order of atoms match")
                else:
                    print("Order of atoms does not match")
            except:
                same_order = False
                print("Could not compare order of atoms, assuming they do not match")
                
            # If the order of the atoms does not match, we try to compute the RMSD using the atom names
            if not same_order:
                try:
                    # Create the RMSD object
                    rmsd_calculator = rms.RMSD(
                        prep.select_atoms(f"resid {resi.resid} and protein and not resname ACE NME NMA"),  # Reference group
                        prev.select_atoms(f"resid {resj.resid} and protein and not resname ACE NME NMA"),  # Trajectory group
                        center=False,  # Do not center coordinates
                        superposition=False  # Disable alignment 
                    )
                    # Compute the RMSD
                    rmsd_calculator.run()
                    if rmsd_calculator.rmsd[0] > 1.0:
                        print("WARNING: RMSD between residues is >1 A")
                except:
                    print("Could not compute RMSD between residues")
                
            if same_order:
                arr = np.abs(resi.atoms.positions - resj.atoms.positions)
                if np.any(arr > 1.0):
                    print("WARNING: >1 A deviation in residue: ", resi.resname, resi.resid)

    prep_nonprot = prep.select_atoms('not protein and not resname NMA')
    prev_nonprot = prev.select_atoms('not protein and not resname NMA')

    # Check that the number of atoms match for non-protein residues
    if len(prep_nonprot.atoms) != len(prev_nonprot.atoms):
        print("WARNING: Number of non protein atoms does not match")

    # Inspect disulfide bridges
    # first guess bonds
    prep.atoms.guess_bonds()
    sgs = prep.select_atoms('name SG and resname CYS CYX')
    distances = distance_array(sgs.atoms.positions, sgs.atoms.positions)

    for i in range(len(sgs)):
        for j in range(len(sgs)):
            # skip self
            if i == j:
                continue
            else:
                if distances[i, j] < 2.5:
                    if sgs.atoms[j] not in sgs.atoms[i].bonded_atoms:
                        if 'H' in sgs.atoms[i].bonded_atoms.elements:
                            print(f"ERROR: possible missed bridge between {sgs.atoms[j]} and {sgs.atoms[i]}")
                        else:
                            print(f"WARNING: unbonded disulfide bridge found: {sgs.atoms[j]} and {sgs.atoms[i]}")
                    else:
                        print(f"LOG: disulfide bridge found {sgs.atoms[j]} and {sgs.atoms[i]}")


if __name__ == "__main__":
    run()
