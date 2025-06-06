import string
import click
import pathlib
import logging
import warnings
import os
import json
from functools import partial
from openff.units import unit
import openfe
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol
from rdkit import Chem
import kartograf
from kartograf.filters import (
    filter_ringbreak_changes,
    filter_ringsize_changes,
    filter_whole_rings_only,
)
from openff.toolkit import (
    RDKitToolkitWrapper, AmberToolsToolkitWrapper
)
from openff.toolkit.utils.toolkit_registry import (
    toolkit_registry_manager, ToolkitRegistry
)
import gufe
from gufe import tokenization


logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", message="Partial charges have been provided, these will preferentially be used instead of generating new partial charges")
amber_rdkit = ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])


def gen_charges(smc):
    """
    Generate AM1BCC partial charges for a SmallMoleculeComponent using
    the input conformer and antechamber as backend.
    """
    print(f"INFO: generating partial charges for ligand {smc.name} -- this may be slow")
    offmol = smc.to_openff()
    with toolkit_registry_manager(amber_rdkit):
        offmol.assign_partial_charges(
            partial_charge_method="am1bcc",
            use_conformers=offmol.conformers
        )
    return openfe.SmallMoleculeComponent.from_openff(offmol)


def gen_ligand_network(smcs, atom_mapper='kartograf', max_3d=0.95):
    """
    Creates the Lomap ligand network using an atom mapper 
    (by default the KartografAtomMapper) and the Lomap scorer.

    Parameters
    ----------
    smcs : list[SmallMoleculeComponents]
      List of SmallMoleculeComponents of the ligands.
    atom_mapper : str
      The atom mapper to use for the ligand network generation.
    max_3d : float
        The maximum distance between 3D coordinates to map atoms.

    Returns
    -------
    openfe.LigandNetwork
      The Lomap generated LigandNetwork.
    """
    print("INFO: Generating Lomap Network")
    mapping_filters = [
        filter_ringbreak_changes,  # default
        filter_ringsize_changes,  # default
        filter_whole_rings_only,  # default
    ]
    
    if atom_mapper == 'kartograf':
        print("INFO: Using Kartograf as the atom mapper")
        mapper = kartograf.KartografAtomMapper(
            atom_max_distance=max_3d,
            atom_map_hydrogens=True,
            additional_mapping_filter_functions=mapping_filters,
        )
    elif atom_mapper == 'lomap':
        print("INFO: Using Lomap as the atom mapper")
        mapper = openfe.LomapAtomMapper(
            max3d = max_3d, 
            element_change = True
        )
    else:
        raise ValueError(f"Unknown atom mapper: {atom_mapper}. Please use 'kartograf' or 'lomap'.")
    
    # TODO: Change this after LOMAP PR gets merged
    # scorer = openfe.lomap_scorers.default_lomap_score
    scorer = partial(openfe.lomap_scorers.default_lomap_score, charge_changes_score=0.1)
    ligand_network = openfe.ligand_network_planning.generate_lomap_network(
        molecules=smcs, mappers=mapper, scorer=scorer)
    # Raise an error if the network is not connected
    if not ligand_network.is_connected():
        errormsg = ('Error encountered when creating the ligand network. The network'
                    ' is unconnected.')
        raise ValueError(errormsg)
    return ligand_network


def get_alchemical_charge_difference(mapping) -> int:
    """
    Checks and returns the difference in formal charge between state A and B.

    Parameters
    ----------
    mapping : dict[str, ComponentMapping]
      Dictionary of mappings between transforming components.

    Returns
    -------
    int
      The formal charge difference between states A and B.
      This is defined as sum(charge state A) - sum(charge state B)
    """
    chg_A = Chem.rdmolops.GetFormalCharge(
        mapping.componentA.to_rdkit()
    )
    chg_B = Chem.rdmolops.GetFormalCharge(
        mapping.componentB.to_rdkit()
    )

    return chg_A - chg_B


def get_settings(min_steps, equil_length, small_molecule_forcefield='openff-2.2.0', hydrogen_mass=3.0):
    """
    Utility method for getting RFEProtocol settings for non charge changing
    transformations.
    These settings mostly follow defaults but use the newest OpenFF 2.2.
    
    Parameters
    ----------
    
    min_steps : int
        Number of minimization steps to perform.
    equil_length : float
        Length of the equilibration in ns.
    small_molecule_forcefield : str
        The name of the small molecule forcefield to use.
    hydrogen_mass : float
        The mass of the hydrogen atoms in the system.
    """
    # Are there additional settings we should specify here?
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Should we use this new OpenFF version or the default?
    settings.forcefield_settings.small_molecule_forcefield = small_molecule_forcefield
    settings.forcefield_settings.hydrogen_mass = hydrogen_mass
    # Set the number of minimization steps
    settings.simulation_settings.minimization_steps = min_steps
    # Set the length of the equilibration
    settings.simulation_settings.equilibration_length = equil_length * unit.nanosecond
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


def get_settings_charge_changes(min_steps, equil_length, small_molecule_forcefield='openff-2.2.0', hydrogen_mass=3.0):
    """
    Utility method for getting RFEProtocol settings for charge changing
    transformations.

    These settings mostly follow defaults but use longer
    simulation times, more lambda windows and an alchemical ion.

    Parameters
    ----------
    
    min_steps : int
        Number of minimization steps to perform.
    equil_length : float
        Length of the equilibration in ns.
    small_molecule_forcefield : str
        The name of the small molecule forcefield to use.
    hydrogen_mass : float
        The mass of the hydrogen atoms in the system.
    """
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Should we use this new OpenFF version or the default?
    settings.forcefield_settings.small_molecule_forcefield = small_molecule_forcefield
    settings.forcefield_settings.hydrogen_mass = hydrogen_mass
    settings.alchemical_settings.explicit_charge_correction = True
    # Set the number of minimization steps
    settings.simulation_settings.minimization_steps = min_steps
    # Set the length of the equilibration
    settings.simulation_settings.equilibration_length = equil_length * unit.nanosecond
    settings.simulation_settings.production_length = 20 * unit.nanosecond
    settings.simulation_settings.n_replicas = 22
    settings.lambda_settings.lambda_windows = 22
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


@click.command
@click.option(
    '--ligands',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared SDF file containing the ligands",
)
@click.option(
    '--pdb',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared PDB file of the protein",
)
@click.option(
    '--cofactors',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=None,
    help="Path to the prepared cofactors SDF file (optional)",
)
@click.option(
    '--min_steps',
    type=int,
    default=5000,
    help="Number of minimization steps to perform",
)
@click.option(
    '--equil_length',
    type=float,
    default=1.0,
    help="Length of the equilibration in ns",
)
@click.option(
    '--ligands_ff',
    type=str,
    default='openff-2.2.0',
    help="The name of the small molecule forcefield to use",
)
@click.option(
    '--hydrogen_mass',
    type=float,
    default=3.0,
    help="The mass of the hydrogen atoms in the system",
)
@click.option(
    '--atom_mapper', 
    type=click.Choice(['lomap', 'kartograf']),
    default='kartograf',
    help="The atom mapper to use for the ligand network generation",
)
@click.option(
    '--max_3d',
    type=float,
    default=0.95,
    help="The maximum distance between 3D coordinates to map atoms",
)
@click.option(
    '--output',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path('alchemicalNetwork'),
    help="Directory name in which to store the transformation json files",
)
def run_inputs(ligands, pdb, cofactors, min_steps, equil_length, ligands_ff, hydrogen_mass, atom_mapper, max_3d, output):
    """
    Generate run json files for RBFE calculations

    Parameters
    ----------
    ligands : pathlib.Path
      A Path to a ligands SDF.
    pdb : pathlib.Path
      A Path to a protein PDB file.
    cofactors : Optional[pathlib.Path]
      A Path to an SDF file containing the system's cofactors.
    min_steps : int
        Number of minimization steps to perform.
    equil_length : float
        Length of the equilibration in ns.
    ligands_ff : str
        The name of the small molecule forcefield to use.
    hydrogen_mass : float
        The mass of the hydrogen atoms in the system.
    atom_mapper : str
        The atom mapper to use for the ligand network generation.
    output: pathlib.Path
      A Path to a directory where the transformation json files
      and ligand network graphml file will be stored into.
    """
    # Create the output directory -- default to alchemicalNetwork, fail if it exists
    output.mkdir(exist_ok=False, parents=True)
    
    # Create the small molecule components of the ligands
    rdmols = [mol for mol in Chem.SDMolSupplier(str(ligands), removeHs=False)]
    smcs = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in rdmols]
    # Generate the partial charges
    logger.info("Generating partial charges for ligands")
    smcs = [gen_charges(smc) for smc in smcs]

    # Create ligand network
    ligand_network = gen_ligand_network(smcs, atom_mapper=atom_mapper, max_3d=max_3d)

    # Store the ligand network as a graphml file
    with open(output / "ligand_network.graphml", mode='w') as f:
        f.write(ligand_network.to_graphml())

    # Create the solvent and protein components
    solv = openfe.SolventComponent()
    prot = openfe.ProteinComponent.from_pdb_file(str(pdb))

    # If we have cofactors, load them in and assign partial charges
    if cofactors is not None:
        cofactors_smc = [gen_charges(openfe.SmallMoleculeComponent(m))
                         for m in Chem.SDMolSupplier(str(cofactors), removeHs=False)]

    # Create the AlchemicalTransformations, and storing them to an AlchemicalNetwork
    transformations = []
    for mapping in ligand_network.edges:
        # Get different settings depending on whether the transformation
        # involves a change in net charge
        charge_difference = get_alchemical_charge_difference(mapping)
        if abs(charge_difference) > 1e-3:
            # Raise a warning that a charge changing transformation is included
            # in the network
            wmsg = ("Charge changing transformation between ligands "
                    f"{mapping.componentA.name} and {mapping.componentB.name}. "
                    "A more expensive protocol with 22 lambda windows, sampled "
                    "for 20 ns each, will be used here.")
            warnings.warn(wmsg)
            # Get settings for charge changing transformations
            rfe_settings = get_settings_charge_changes(min_steps, equil_length, ligands_ff, hydrogen_mass)
        else:
            rfe_settings = get_settings(min_steps, equil_length, ligands_ff, hydrogen_mass)
            
        # Print settings into a different file for this edge
        with open(output / f"{mapping.componentA.name}_{mapping.componentB.name}_settings.txt", "w") as f:
            # Write all the settings
            for key, value in rfe_settings.dict().items():
                f.write(f"{key}: {value}\n")
                
        for leg in ['solvent', 'complex']:
            # use the solvent and protein created above
            sysA_dict = {'ligand': mapping.componentA,
                         'solvent': solv}
            sysB_dict = {'ligand': mapping.componentB,
                         'solvent': solv}

            if leg == 'complex':
                sysA_dict['protein'] = prot
                sysB_dict['protein'] = prot
                if cofactors is not None:

                    for cofactor, entry in zip(cofactors_smc,
                                               string.ascii_lowercase):
                        cofactor_name = f"cofactor_{entry}"
                        sysA_dict[cofactor_name] = cofactor
                        sysB_dict[cofactor_name] = cofactor

            sysA = openfe.ChemicalSystem(sysA_dict)
            sysB = openfe.ChemicalSystem(sysB_dict)

            name = (f"{leg}_{mapping.componentA.name}_"
                    f"{mapping.componentB.name}")

            rbfe_protocol = RelativeHybridTopologyProtocol(settings=rfe_settings)
            transformation = openfe.Transformation(
                stateA=sysA,
                stateB=sysB,
                mapping=mapping,
                protocol=rbfe_protocol,
                name=name
            )
            transformations.append(transformation)

    # Create the alchemical network and write it to disk
    alchemical_network = openfe.AlchemicalNetwork(transformations)
    alchemical_network_json_fp = output / "alchemical_network.json"
    json.dump(
        alchemical_network.to_dict(),
        alchemical_network_json_fp.open(mode="w"),
        cls=tokenization.JSON_HANDLER.encoder
    )

    # Write out each transformation
    # Create a subdirectory for the transformations
    transforms_dir = pathlib.Path(output / "transformations")
    transforms_dir.mkdir(exist_ok=True, parents=True)

    for transform in alchemical_network.edges:
        transform.dump(transforms_dir / f"{transform.name}.json")


if __name__ == "__main__":
    run_inputs()
