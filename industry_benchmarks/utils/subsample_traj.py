import MDAnalysis as mda
import numpy as np
import pathlib
import argparse

from openfe_analysis import FEReader

def subsample_traj(simulation, hybrid_system_pdb, lambda_windows, outfile):
    
    for i in range(0, lambda_windows):
        
        u = mda.Universe(hybrid_system_pdb, simulation, format=FEReader, state_id=i)
        
        # Show details of the trajectory
        print(f"Number of frames: {len(u.trajectory)}")
        #print(f"Number of atoms: {len(u.atoms)}")
        
        out_traj = pathlib.Path(f"{outfile}_{i}.xtc")
        with mda.Writer(str(out_traj), n_atoms=len(u.atoms)) as w:
            for ts in u.trajectory:
                w.write(u.atoms)
                
                
def main():
    parser = argparse.ArgumentParser(
        description="Subsample a trajectory for each lambda window"
    )
    parser.add_argument(
        "--traj_path", type=str, help="Path to the trajectory file"
    )
    parser.add_argument(
        "--top_path", type=str, help="Path to the topology file"
    )
    parser.add_argument(
        "--lambda_windows", type=int, help="Number of lambda windows"
    )
    parser.add_argument(
        "--outfile", type=str, help="Output file name"
    )

    args = parser.parse_args()
    
    subsample_traj(args.traj_path, args.top_path, args.lambda_windows, args.outfile)

if __name__ == "__main__":
    main()