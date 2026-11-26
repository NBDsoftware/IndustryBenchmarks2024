# Running a new system

How to run an RBFE benchmark on a new system with the scripts in this folder, following the
[IB2024 protocol](https://industrybenchmarks2024.readthedocs.io/en/latest/).

Example SLURM scripts for every step are in `/shared/work/pnavarro/projects/2024_ALMIRALL/Benchmarks`
(see `CRBN/` and `JAK3/`). Copy them and change the paths.

## Before you start

- **Environment** -> a conda environment with OpenFE 1.4 (see `environment.yml`). In the SLURM scripts: `module load Miniconda3` and `source activate /path/to/openfe/env`.
- **Inputs** -> a prepared protein PDB, a ligands SDF with 3D poses in the binding site (unique,
  anonymised names), and optionally a cofactors SDF.
- **Suggested layout**, one folder per system:

  ```
  SYSTEM/
  ├── input/                          # protein.pdb, ligands.sdf, cofactors.sdf
  ├── rbfe_network_planning_lomap/    # planning job
  ├── rbfe_network_planning_kartograf/
  ├── transformations_lomap/          # planning output
  ├── transformations_kartograf/
  ├── campaign_results_0/             # one folder per repeat
  ├── campaign_results_1/
  └── campaign_results_2/
  ```

## Workflow

1. [Validate the inputs](#1-validate-the-inputs).
2. [Align the ligands](#2-align-the-ligands) to a reference and check the result. Mantaining or not the 3D conformations.
3. [Plan the networks](#3-plan-the-networks) with Lomap and Kartograf, using the best aligned series.
4. [Compare and choose a network](#4-compare-and-choose-a-network).
5. [Run the campaign](#5-run-the-campaign).
6. [Analyse the results](#6-analyse-the-results).

## 1. Validate the inputs

Run a short MD of the protein (and cofactors) to catch problems such as bad residue or terminal
names. Needs a GPU.

```bash
python input_validation.py --pdb protein.pdb [--cofactors cofactors.sdf]
```

- Success: the output ends with `SIMULATION COMPLETE`.
- Optional: `maint/check_protein_prep.py --prepared prepared.pdb --original original.pdb` compares
  the prepared protein against the original structure.

## 2. Align the ligands

Both atom mappers only map atoms that are closer than `--max_3d` (0.95 A by default) in 3D. Poorly
overlapping poses give larger transforming regions, so the alignment matters, mostly for Kartograf.

Use [notebooks/align_series.ipynb](notebooks/align_series.ipynb):

1. Set `ligands_path`, `ref_ligand_path` (crystal ligand or best docked pose), and `output_path`.
2. Run the alignment. It finds the MCS of the whole series plus the reference, and rigidly aligns each
   ligand onto the reference over the MCS atoms. Conformations are not changed.
3. Check the evaluation section:
   - Pairwise MCS RMSD (in place, no fitting) before and after the alignment. Lower is better.
   - Heatmaps and histogram of the pairwise RMSD.
   - Mean RMSD per ligand, to spot outliers.
4. Look at the aligned poses in the binding site. Rigid alignment can create clashes with the
   protein.

Keep the series (original or aligned) with the best overlap and no clashes. For example, in JAK3 the
re-aligned poses clashed with the protein, so the docked poses were used.

## 3. Plan the networks

Plan one network per atom mapper with the chosen series, so you can compare them in the next step.
Runs on CPU (8 tasks, under 2 h in our systems).

```bash
python plan_rbfe_network.py \
    --pdb protein.pdb \
    --ligands ligands.sdf \
    [--cofactors cofactors.sdf] \
    --atom_mapper lomap \
    --output ../transformations_lomap
```

What it does:

- Assigns AM1BCC partial charges to the ligands once, so all repeats share them.
- Builds atom mappings for all ligand pairs, scores them, and keeps a minimal connected network
  (Lomap network planner).
- Writes to `--output`:
  - `transformations/` -> one JSON per edge and leg (`complex` and `solvent`), each runs 1 repeat.
  - `ligand_network.graphml` -> the ligand network with its atom mappings.
  - `alchemical_network.json` -> the full alchemical network.

Useful options:

| Option | Default | Notes |
|---|---|---|
| `--atom_mapper` | `kartograf` | `lomap` or `kartograf` |
| `--max_3d` | `0.95` | Max distance (A) between mapped atoms. JAK3 used 3 with Lomap. |
| `--ligands_ff` | `openff-2.2.0` | Small molecule force field |
| `--equil_length` | `1.0` | Equilibration length (ns) |

## 4. Compare and choose a network

Use [notebooks/visualize_network.ipynb](notebooks/visualize_network.ipynb):

1. Set `ligand_network_path` and `ligand_network_path_2` to the two `ligand_network.graphml` files.
2. Draw the network, and a single edge mapping by setting `A_name` and `B_name`.
3. Compare the transforming region of each edge:
   - **red** -> atoms created or destroyed.
   - **blue** -> mapped atoms that change element.
   - Summary per network (mean, median, and max), histograms, and the shared edges.

How to choose:

- Smaller transforming regions tend to converge better, so prefer the network with fewer large edges.
- Sometimes a larger transforming region is needed, to allow groups to rotate and to get enough phase
  space overlap between ligands. Check the largest edges by hand before choosing.
- For example, JAK3 used Lomap because Kartograf gave large transforming regions due to binding mode
  changes.

## 5. Run the campaign

Each JSON is one GPU job: `openfe quickrun` runs the leg, then `results_cleanup.py` shrinks the
result JSON. The total number of jobs is edges x 2 legs x 3 repeats.

1. **Test one edge first.** Run one complex leg (the slowest) to check it works and to estimate the
   time per job against the SLURM time limit.
2. **Launch the repeats.** Copy `run_campaing.sh` and `run_transformation.sl` from the example
   folder, set `transformations_path` and `system_name`, and run it once inside each of
   `campaign_results_0`, `_1`, and `_2`:

   ```bash
   cd campaign_results_0
   bash ../run_campaing.sh
   ```

   Results go to `results/`, working files to `simulations/`.
3. **Handle failed edges** as described in the
   [IB2024 docs](https://industrybenchmarks2024.readthedocs.io/en/latest/public/overview.html#handling-failed-edges):
   rerun failed jobs up to 3 times, drop edges that always fail if they are redundant.
   For NaN errors, `find_max_force.py` lists the atoms with the largest forces in the dumped state.

## 6. Analyse the results

Gather the dG and ddG estimates of all repeats:

```bash
python extras/extract_results.py \
    --results_path campaign_results_0/results \
    --results_path campaign_results_1/results \
    --results_path campaign_results_2/results \
    --input_ligand_network_file transformations_lomap/ligand_network.graphml \
    -o ddG_results.tsv \
    -o_DG dG_results.tsv
```

- `ddG_results.tsv` -> ddG per edge with its uncertainty.
- `dG_results.tsv` -> dG per ligand (MLE). Needs cycles in the network.

Then use [notebooks/plot_results.ipynb](notebooks/plot_results.ipynb) to read the experimental values
from the ligands SDF and plot the correlations against experiment with cinnabar.
