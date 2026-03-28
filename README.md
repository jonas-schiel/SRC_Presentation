#BYU Student Research Conference Presentation, February 2026
##Clinical Trial Dose-Escalation Models

This repository contains my research presentation, as well as the R code used for simulating and comparing several Phase 1 dose-escalation trial designs, including the 3+3, mTPI, mTPI-2, and CRM (Continual Reassessment Method) approaches. The project was developed as part of research with the BYU Statistics Department.

---

## Dependencies

- R (>= 4.0)
- `rstan` or `cmdstanr` (for CRM model)
- `here` (for portable file paths)
- Any additional packages are loaded at the top of each script

Install the `here` package if you don't have it:
```r
install.packages("here")
```

---

## Repository Structure

```
/
├── README.md
├── presentation.pdf          # Presentation slides (PDF)
├── presentation.key          # Presentation slides (Keynote)
├── plots/                    # Scripts for generating plots used in presentation
└── trial_design/             # Scripts for simulating each trial design once, and simulating trials at scale
```

### `plots/`

Contains scripts for visualizing simulation results. Run these after generating simulation data from the scripts in `trial_design/`.

| File | Description |
|---|---|
| `LargeSimPlots.R` | Generates plots for large-scale simulations of the 3+3, mTPI, and mTPI-2 trial designs |
| `SingleTrialPlots.R` | Generates plots for individual trial runs of the 3+3, mTPI, and mTPI-2 trial designs |
| `crmplots.R` | Generates plots specific to the CRM model |

### `trial_design/`

Contains all simulation logic. Organized into subfolders by design type.

| File / Folder | Description |
|---|---|
| `LargeTrialSim.R` | Runs large-scale simulations for the 3+3, mTPI, and mTPI-2 trial designs |
| `base_trial_sims/` | Individual simulation functions for the 3+3, mTPI, and mTPI-2 trial designs |
| `crm_design/` | CRM model implementation using Stan |

#### `trial_design/base_trial_sims/`

| File | Description |
|---|---|
| `3by3_function.R` | Simulates a single 3+3 design trial |
| `mtpi_function.R` | Simulates a single mTPI design trial |
| `mtpi2_function.R` | Simulates a single mTPI-2 design trial |

#### `trial_design/crm_design/`

Contains the Stan model and associated R code for the Continual Reassessment Method (CRM). The Stan model is compiled at runtime via `rstan` or `cmdstanr`.

---

## Usage

1. Clone the repository and open the `.Rproj` file in RStudio — this ensures `here()` resolves paths correctly from the repo root.
2. To simulate individual trials, source the relevant function file from `base_trial_sims/` or run the CRM scripts in `crm_design/`.
3. To run large-scale simulations, run `LargeTrialSim.R` from the `trial_design/` folder.
4. To generate plots, run the relevant script from the `plots/` folder after simulation data has been generated.

---

## Notes

- All file paths use the `here` package for portability. As long as you open the project via the `.Rproj` file, paths should resolve correctly on any machine.
- The CRM model requires a working Stan installation. See the [RStan getting started guide](https://mc-stan.org/rstan/) for setup instructions.
