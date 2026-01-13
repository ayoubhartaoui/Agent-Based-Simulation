# Master Thesis – Simulation and Analysis Code

This repository contains the simulation and analysis code developed for my Master thesis. The project implements an evolutionary simulation model of social interactions (conflict, trade, and avoidance) under different environmental variability regimes, together with the analysis pipeline used to produce the results presented in the thesis.

Large simulation outputs are not included in this repository due to their size.

---

## Repository structure

```
.
├─ Simulation/
│  ├─ run_simulation.jl
│  ├─ Project.toml
│  └─ Manifest.toml
│
├─ Analysis/
│  └─ analysis.R
│
├─ output/
│  └─ .gitkeep
│
├─ README.md
├─ .gitignore
└─ LICENSE
```

---

## Installation (Julia)

Instantiate the Julia environment used for the simulations:

```bash
julia --project=Simulation -e 'using Pkg; Pkg.instantiate()'
```

---

## Running the simulations

The main simulation script (covering all environmental regimes and trade rules studied in the thesis) is:

```
Simulation/run_simulation.jl
```

Run the simulations with:

```bash
julia --project=Simulation Simulation/run_simulation.jl
```

All output files are written locally and are not tracked by git.

---

## ⚠️ Output size and storage considerations (READ THIS)

Each full simulation run produces extremely large output files.

When run with the default parameters used in the Master thesis:

- Total output size can reach **~20 GB per run**
- Outputs consist of multiple CSV files containing:
  - full time series
  - replicate-level summaries
  - within-population samples over the final time window

These files are intentionally not included in the repository and are excluded via `.gitignore`.

### Strong recommendations

- Ensure you have sufficient disk space before running the simulations.
- Prefer running the code:
  - on a machine with large local storage, or
  - on a server with adequate user storage allocation.

During the thesis, simulations were stored using the **University of Lausanne OneDrive allocation (100 GB)**.

### HPC usage warning

If you plan to run this code on an HPC cluster:

- Be mindful of storage quotas.
- Avoid writing large CSV files to scratch or quota-limited directories.
- Consider:
  - reducing the number of generations or replicates for test runs,
  - redirecting outputs to long-term storage,
  - or running simulations selectively rather than in bulk.

Blindly running the default configuration on an HPC system may exceed storage limits.

---

## Analysis

After simulations have completed and CSV files are available locally, the analysis can be run with:

```bash
Rscript Analysis/analysis.R
```

---

## License

See the `LICENSE` file for details.

