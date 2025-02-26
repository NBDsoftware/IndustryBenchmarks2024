import argparse
import pathlib
from rdkit import Chem
from openff.units import unit
import gufe
import openfe
from openfe.protocols.openmm_md.plain_md_methods import PlainMDProtocol
from openfe import ChemicalSystem, ProteinComponent, SmallMoleculeComponent, SolventComponent

def parse_args():
    parser = argparse.ArgumentParser(description="Run an MD simulation with a given protein and ligand.")
    parser.add_argument("--pdb", required=True, help="Path to the protein PDB file.")
    parser.add_argument("--ligand", required=False, help="Path to the ligand SDF file.")
    parser.add_argument("--ligand_name", required=False, help="Name of the ligand to use in the simulation.")
    parser.add_argument("--cofactors", required=False, help="Path to the cofactors SDF file.")
    parser.add_argument("--output", required=False, help="Path to the output directory.")
    return parser.parse_args()

def get_settings():

    settings = PlainMDProtocol.default_settings()
    settings.simulation_settings.equilibration_length_nvt = 1 * unit.nanosecond
    settings.simulation_settings.equilibration_length = 1 * unit.nanosecond
    settings.simulation_settings.production_length = 5 * unit.nanosecond
    
    return settings
    
def load_ligand(ligand_path, ligand_name):
    ligands_sdf = Chem.SDMolSupplier(ligand_path, removeHs=False)
    ligands = [SmallMoleculeComponent(mol) for mol in ligands_sdf if mol is not None]
    
    if not ligands:
        raise ValueError(f"No valid molecules found in {ligand_path}.")
    
    ligand = next((lig for lig in ligands if lig.name == ligand_name), None)
    if ligand is None:
        raise ValueError(f"Ligand '{ligand_name}' not found in {ligand_path}.")
    
    return ligand

def load_cofactors(cofactors_path):
    cofactors_sdf = Chem.SDMolSupplier(cofactors_path, removeHs=False)
    cofactors = [SmallMoleculeComponent(mol) for mol in cofactors_sdf if mol is not None]
    
    if not cofactors:
        raise ValueError(f"No valid molecules found in {cofactors_path}.")
    
    return cofactors

def run_md(dag, output):
    """
    Run a DAG and check it was ok.

    Parameters
    ----------
    dag : openfe.ProtocolDAG
      A ProtocolDAG to execute.
    output : str or None
        The directory to run the simulation in.

    Raises
    ------
    AssertionError
      If any of the simulation Units failed.
    """
    # Add output directory if any
    if output:
        output = pathlib.Path(output)
    else:
        output = pathlib.Path.cwd() / "output"
        
    # Create the output directory if it does not exist
    output.mkdir(parents=True, exist_ok=True)
        
    # Run the DAG
    dagres = gufe.protocols.execute_DAG(
        dag,
        shared_basedir=output,
        scratch_basedir=output,
        keep_shared=True,
        raise_error=True,
        n_retries=0,
    )

    # If everything is good then tell the user
    assert dagres.ok()
    print("SIMULATION COMPLETE")
    
def main():
    args = parse_args()
    
    # Create components dictionary
    components_dict = {}
    
    # Load protein
    components_dict['protein'] = ProteinComponent.from_pdb_file(args.pdb, name='protein')
    name = components_dict['protein'].name
    
    # Load solvent
    components_dict['solvent'] = SolventComponent(ion_concentration=0.15 * unit.molar)
    
    # Load ligand if any
    if args.ligand:
        ligand = load_ligand(args.ligand, args.ligand_name)
        components_dict['ligand'] = ligand
        name = f"{name}_{ligand.name}"
    
    # Load cofactors if any
    if args.cofactors:
        cofactors = load_cofactors(args.cofactors)
        components_dict['cofactors'] = cofactors
        
    # Create the ChemicalSystem
    system = ChemicalSystem(components_dict, name=name)
    
    # Create and execute MD simulation
    settings = get_settings()
    protocol = PlainMDProtocol(settings=settings)
    dag = protocol.create(stateA=system, stateB=system, mapping=None)
        
    run_md(dag, args.output)
    
    print("MD simulation completed successfully.")

if __name__ == "__main__":
    main()
