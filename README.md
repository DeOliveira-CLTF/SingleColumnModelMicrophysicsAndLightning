# Evaluating lightning parameterizations in a single-column model

Code and reproducibility materials for the manuscript:

> **Evaluating lightning parameterizations and their sensitivity to cloud
> microphysics in a single-column model: a severe-convection case study**

The study evaluates lightning-rate parameterizations driven by simulated
cloud and precipitation properties. The experiments use a one-dimensional
single-column model (SCM), two cloud microphysics schemes (Thompson and
WSM6), and three lightning parameterizations. The main scientific question is
how the lightning predictions and their sensitivity to microphysics change
when the same atmospheric case is simulated with different microphysics.

The repository contains the model-side diagnostics, observation processing,
calibration, ensemble experiments, evaluation, and publication figures. It
does **not** redistribute restricted lightning observations, private data, or
credentials.

## Study design

Two independent cases are used:

| Case | Role | Period | Use |
|---|---|---|---|
| Argentina 2019 | Development/calibration | 28–29 October 2019 | Fit the lightning coefficients using deterministic controls |
| Argentina 2018 | Independent evaluation | 13–14 December 2018 | Apply the frozen coefficients without recalibration |

The 2018 observations are never used to fit the coefficients. The active
evaluation uses one-minute lightning observations and reports lightning rates
in flashes min⁻¹. The model ensemble contains 40 perturbed members plus one
control for each microphysics/coupling configuration.

The principal configurations are:

- Thompson and WSM6 cloud microphysics;
- C0, CL, and CF MYNN coupling options;
- CAPE × precipitation (`MUCAPE × P`);
- an ice–graupel product (`I × G`);
- the adapted McCaul F3 formulation using prescribed large-scale vertical
  motion.

The physical SCM, forcing, and observation data are kept separate from the
plotting and statistical code. Radiosondes are regional reference profiles;
they are not used as SCM initial or boundary conditions and are not assumed
to be collocated with the fixed model column.

## Repository map

| Path | Purpose |
|---|---|
| `src/` | Fortran SCM and model diagnostics |
| `python/scm_lightning/` | Python package for forcing, observations, lightning diagnostics, ensembles, calibration, evaluation, and plotting |
| `scripts/` | Public command-line entry points |
| `configs/cases/` | Case definitions, dates, domains, and data paths |
| `configs/evaluation/` | Calibration and validation protocols |
| `configs/experiments/` | Sensitivity and experiment designs |
| `data/raw/` | Local input-data cache; normally not distributed |
| `runs/` | Model executions and member outputs |
| `results/` | Tables, diagnostics, manifests, and publication figures |
| `docs/` | Data sources, methodology, installation, and reproducibility notes |
| `tests/` | Automated tests and scientific contract checks |
| `external/DEPHY-SCM/` | DEPHY-SCM dependency/submodule |
| `build/` | Local compiled executables; ignored by Git |

Important Python subpackages:

- `forcing/`: ERA5 loading, normalization, spatial subsetting, and forcing
  contracts;
- `observations/`: lightning-event and radiosonde processing;
- `lightning/`: lightning-rate predictors and physical diagnostics;
- `ensemble/`: controls, perturbations, member preparation, and execution;
- `calibration/`: coefficient fitting and calibration metadata;
- `evaluation/`: independent evaluation, metrics, audits, and tables;
- `plotting/`: article and supplementary figures.

## Requirements

- Python 3.10 or newer;
- Python 3.12 is the tested Conda environment;
- Fortran compiler and NetCDF-Fortran for SCM compilation;
- Linux or WSL is recommended for the Fortran build;
- access to the required ERA5 and project-observation files.

The Python dependencies are declared in [`pyproject.toml`](pyproject.toml) and
[`environment.yml`](environment.yml). Create the Conda environment with:

```bash
conda env create -f environment.yml
conda activate scm-lightning
python -m pip install -e ".[test]"
```

Alternatively, with an existing Python environment:

```bash
python -m pip install -e ".[test]"
```

Do not commit CDS, ENTLN, Google Drive, or other private credentials. Configure
access outside the repository and keep private input data in the paths
referenced by the case YAML files.

## Quick start

Run every command from the repository root. On Windows PowerShell, replace
`python` with the selected interpreter if necessary. On Linux/WSL, use the
same commands from the mounted repository.

Inspect the available stages before running a workflow:

```bash
python scripts/reproduce.py --help
python scripts/download.py --help
python scripts/prepare_forcing.py --help
python scripts/run_scm.py --help
python scripts/process_lightning.py --help
python scripts/calibrate.py --help
python scripts/evaluate.py --help
python scripts/plot.py --help
```

Run the automated tests:

```bash
python -m pytest tests -q
```

Check the current reproducibility state:

```bash
python scripts/reproduce.py status
```

## Reproduction workflow

The stages are intentionally separate so that data preparation, model
execution, calibration, evaluation, and plotting can be audited independently.
Do not skip the provenance/manifests when creating a new result tree.

### 1. Prepare data and forcing

```bash
python scripts/download.py --help
python scripts/download.py boundaries
python scripts/prepare_forcing.py --help
```

The exact data locations are defined in the case configuration. Review
[`docs/data_sources.md`](docs/data_sources.md) before downloading or copying
large input files.

### 2. Build and run the SCM

Builds must be performed in Linux/WSL with `gfortran` and NetCDF-Fortran:

```bash
make -j4 BUILD=build/bulk_memory_final
```

The public execution stages are:

```bash
python scripts/run_scm.py controls
python scripts/run_scm.py prepare-ensemble
python scripts/run_scm.py ensemble
```

Use the command help and the case configuration before changing an output
directory. A new executable or changed input must use a new output tree; do
not overwrite an audited experiment.

### 3. Process observations and model diagnostics

```bash
python scripts/process_lightning.py profiles
python scripts/process_lightning.py observations
python scripts/process_lightning.py aggregate
python scripts/process_lightning.py diagnostics
```

The observed one-minute lightning products are generated from the event
inventory, not by dividing 15-minute counts. The resulting rates are in
flashes min⁻¹; cumulative products remain counts of flashes.

### 4. Calibrate only with the 2019 case

```bash
python scripts/calibrate.py independent
python scripts/calibrate.py fit
python scripts/calibrate.py adaptations
```

Calibration coefficients are stored with their manifests under
`results/argentina_20191028/`. The 2018 case must not be used for fitting.
See [`docs/independent_events.md`](docs/independent_events.md) for the active
calibration/evaluation protocol.

### 5. Evaluate the independent 2018 case

```bash
python scripts/evaluate.py independent
python scripts/evaluate.py article
python scripts/evaluate.py contracts
python scripts/evaluate.py migration
```

The evaluation stages produce metrics, contract checks, provenance manifests,
and tables under `results/argentina_20181213/`.

### 6. Generate article figures

```bash
python scripts/plot.py article
python scripts/plot.py sensitivity
python scripts/plot.py all
```

Publication products are written to:

```text
results/argentina_20181213/article_bulk_memory_final/
```

The figure exporter writes PNG (300 dpi), PDF, and SVG versions where
applicable. The principal products include:

- `figures/Figure_01_case_environment.*`: ERA5 vertical motion, GOES-16
  infrared, and regional radiosonde profiles;
- `figures/Figure_02_observed_lightning.*`: one-minute ENTLN lightning rate
  and cumulative activity;
- `figures/Figure_03_ensemble_predictions.*`: ensemble lightning predictions;
- `figures/Figure_04_performance_sensitivity.*`: performance metrics;
- `figures/Figure_05_microphysical_response.*`: hydrometeor and precipitation
  responses;
- `figures/Figure_06_*`: precipitation comparison;
- `supplementary/S03_*`, `S04_*`, `S_CAPE_*`, and `S2019_*`: supplementary
  analyses and the 2019 environmental characterization.

The exact filenames are recorded in the result manifests. Do not interpret
the presence of an old file in `results/` as evidence that it belongs to the
active protocol; historical products are kept in `legacy/` or archived trees.

## Deterministic vertical-motion sensitivity

The updraft experiment changes only the prescribed large-scale vertical-motion
perturbation while retaining the fixed meteorological design:

```bash
python scripts/run_scm.py sensitivity design
python scripts/run_scm.py sensitivity run
python scripts/process_lightning.py sensitivity
```

Build its executable separately:

```bash
make -j4 BUILD=build/updraft_sensitivity
```

The 18 deterministic cases are stored under
`runs/argentina_20181213/updraft_sensitivity/`. This experiment is not
multiplied by the 40-member initial-condition ensemble and is not calibrated
against lightning observations.

## Scientific conventions

- Native ERA5 omega is retained in Pa s⁻¹; negative omega denotes ascent.
- The SCM diagnostic relation is `w = -omega / (rho g)`.
- `mf_w` is plume mass flux multiplied by mean plume velocity at interfaces.
- Inactive EDMF plumes are represented by `w_EDMF = NaN`; plume velocities are
  not summed.
- The main physics uses `mynn_edmf_type=0` and `mynn_uedmf_mp=0`.
- `rainnc` is used for total accumulated precipitation; `snownc` is not added.
- CAPE × P uses the simulated-profile MetPy pseudoadiabatic CAPE diagnostic;
  it is not presented as the isentropic CAPE formulation from Romps.
- Adapted one-dimensional predictors are identified explicitly in tables and
  metadata; they are not claimed to be identical to their published
  three-dimensional formulations.

For the detailed definitions and limitations, see
[`docs/methodology.md`](docs/methodology.md) and
[`docs/lightning_schemes.md`](docs/lightning_schemes.md).

## Reproducibility and data policy

Every audited result tree should contain manifests with source paths, hashes,
configuration identity, and protocol information. Native NetCDF files and
restricted lightning data must remain unchanged and must not be uploaded to
public services.

Useful checks include:

```bash
python scripts/reproduce.py audit
python scripts/evaluate.py contracts
python scripts/evaluate.py migration
```

The repository preserves older analyses for comparison, but those archives are
not interchangeable with the active 2019-calibration/2018-evaluation
workflow. Read [`docs/reproducibility.md`](docs/reproducibility.md) before
combining products from different result trees.

## Related documentation

- [Installation](docs/installation.md)
- [Data sources](docs/data_sources.md)
- [Methodology](docs/methodology.md)
- [Lightning schemes](docs/lightning_schemes.md)
- [Independent-event protocol](docs/independent_events.md)
- [Reproducibility and contracts](docs/reproducibility.md)

## License and citation

The software license is provided in [`LICENSE`](LICENSE). If you use this
repository, cite the associated manuscript:

> *Evaluating lightning parameterizations and their sensitivity to cloud
> microphysics in a single-column model: a severe-convection case study.*

Please also cite the underlying SCM, ERA5, lightning-observation, satellite,
and radiosonde data products according to their respective terms and
documentation.
