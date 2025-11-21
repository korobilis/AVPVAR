# Replication Package: Learning from Crises - A New Class of Time-Varying Parameter VARs with Observable Adaptation

This repository contains the complete replication package for Hardy and Korobilis (2025): "Learning from crises: A new class of time-varying parameter VARs with observable adaptation".

## Overview

This replication package provides all code and data necessary to reproduce the empirical results, including:
- Out-of-sample forecasting exercises for US monthly, Euro Area quarterly, and FRED-QD data
- In-sample analysis and parameter dynamics
- Monte Carlo simulation studies
- All tables and figures from the main paper and online supplement

## Paper Abstract

[Add your paper abstract here]

## Software Requirements

- **MATLAB**: R2024b or R2025a (tested versions)
- **No proprietary MATLAB toolboxes required**
- **Operating System**: Windows, macOS, or Linux
- **Hardware**: Full replication may take several hours depending on your system

## Repository Structure

```
Replication_Package/
│
├── Forecasting/                    # Out-of-sample forecasting applications
│   ├── EuroArea/                   # Euro Area quarterly forecasting
│   │   ├── data/                   # EA data files
│   │   ├── Forecasting_Results/    # Results and figure generation scripts
│   │   ├── functions/              # Model estimation functions
│   │   ├── FREDQD/                 # FRED-QD quarterly application
│   │   └── matlabtoolbox/          # Helper functions
│   │
│   └── US/                         # US monthly forecasting
│       ├── data/                   # US data files
│       ├── Forecasting_Results/    # Results and figure generation scripts
│       ├── functions/              # Model estimation functions
│       └── matlabtoolbox/          # Helper functions
│
├── InSample/                       # In-sample analysis
│   ├── EA/                         # Euro Area in-sample analysis
│   │   ├── data/
│   │   ├── functions/
│   │   └── Insample_results/
│   │
│   └── US/                         # US in-sample analysis
│       ├── data/
│       ├── functions/
│       └── Insample_Results/
│
└── Monte_Carlo/                    # Monte Carlo simulations
    ├── Figures/                    # Generated figures
    ├── functions/                  # DGP and helper functions
    └── MODELS/                     # Model implementations
```

## Quick Start Guide

### 1. Clone the Repository
```bash
git clone https://github.com/[your-username]/[repo-name].git
cd [repo-name]
```

### 2. Main Replication Scripts

Each application has its own main script. Run these from MATLAB in the respective directories:

#### US Monthly Forecasting
```matlab
cd Forecasting/US/
main_forecasting_US
```

#### Euro Area Quarterly Forecasting
```matlab
cd Forecasting/EuroArea/
main_forecasting_EA
```

#### FRED-QD Forecasting
```matlab
cd Forecasting/EuroArea/FREDQD/
main_forecasting_FREDQD
```

#### In-Sample Analysis
```matlab
cd InSample/US/
INSAMPLE_models

cd InSample/EA/
main_INSAMPLE_EA
```

#### Monte Carlo Simulations
```matlab
cd Monte_Carlo/
MONTE_CARLO
```

## Models Implemented

The package includes nine forecasting models:

1. **AVP-VAR**: Adaptive-parameter VAR with observable drivers
2. **CP-VAR**: Constant-parameter VAR
3. **TVP-VAR-EB**: TVP-VAR with empirical Bayes
4. **CP-VAR-SV**: Constant parameters with stochastic volatility
5. **OLS VAR**: Benchmark homoskedastic VAR
6. **VAR-SVO-t**: VAR with outliers and Student-t errors
7. **FAVAR**: Factor-augmented VAR
8. **FAVAR-SV**: FAVAR with stochastic volatility
9. **TVP-VAR-FB**: Full Bayesian TVP-VAR

## Generating Tables and Figures

After running the main scripts, generate paper figures and tables:

### Forecasting Results
```matlab
cd Forecasting/[US or EuroArea]/Forecasting_Results/
Figures_oos                    % Generate MAE and quantile score figures
Appendix_Forecasting_[US/EA]   % Generate MSPE tables
```

### In-Sample Results
```matlab
cd InSample/[US or EA]/
plot_UCSV_TVintercepts         % Intercept comparison figures
plot_residuals                 % Residual diagnostic figures
```

### Monte Carlo Results
```matlab
cd Monte_Carlo/
Figures1_2                     % DGP illustrations
Appendices_Tables_2_3          % Monte Carlo tables
```

## Data Description

### US Monthly Data
- **Variables**: Industrial Production (INDPRO), PCE Price Index (PCEPI), Federal Funds Rate (FEDFUNDS)
- **Sample**: [Add sample period]
- **Transformations**: Log-differences for INDPRO and PCEPI, first difference for FEDFUNDS

### Euro Area Quarterly Data
- **Variables**: Industrial Production (YER), HICP inflation, Short-term rate (STN)
- **Sample**: [Add sample period]
- **Drivers**: Multiple economic indicators including POILU, PCOMU, COMPR, etc.

### FRED-QD Data
- **Variables**: GDP growth (GDPC1), PCE inflation (PCECTPI), Federal Funds Rate
- **Sample**: Quarterly frequency from FRED-QD database
- **Drivers**: S&P 500, permits, spreads, VIX, sentiment indices

## Expected Runtime

- **Single forecasting application (US/EA)**: 2-4 hours
- **FRED-QD application**: 2-3 hours
- **In-sample analysis**: 1-2 hours per dataset
- **Monte Carlo (1000 replications)**: 6-8 hours
- **Total full replication**: ~20 hours

Runtime depends on hardware specifications and MCMC settings.

## Output Files

Results are saved as `.mat` files:
- `oos_[region]_p[1-4].mat`: Out-of-sample forecasting results
- `insample_results_p[1-2].mat`: In-sample estimation results
- `Monte_Carlo_DGP[X]_All_Configurations.mat`: Monte Carlo results

Figures are saved as PDFs in respective `Figures/` or `Appendix_Figures/` folders.

## MCMC Settings

Default settings across applications:
- **Posterior draws saved**: 5,000
- **Burn-in**: 1,000
- **Thinning**: 5
- **Forecast horizons**: h=24 (US monthly), h=8 (EA/FRED-QD quarterly)

## Detailed Documentation

For detailed instructions on each component, please refer to the individual README files:
- `README_Forecasting_US.pdf`: US monthly forecasting details
- `README_Forecasting_EA.pdf`: Euro Area forecasting details
- `README_Forecasting_FREDQD.pdf`: FRED-QD application details
- `README_Insample.pdf`: In-sample analysis details
- `README_Montecarlo.pdf`: Monte Carlo simulation details

## Citation

If you use this code in your research, please cite:

```bibtex
@article{hardy2025learning,
  title={Learning from crises: A new class of time-varying parameter VARs with observable adaptation},
  author={Hardy, Nicolas and Korobilis, Dimitris},
  year={2025},
  journal={[Journal Name]},
  volume={},
  pages={}
}
```

## Authors

- **Nicolas Hardy** - nicolas.hardy.hernandez@gmail.com
- **Dimitris Korobilis** - dikorobilis@googlemail.com

## License

[Add license information - e.g., MIT, GPL, etc.]

## Contact

For questions or issues with the replication package, please contact:
- Nicolas Hardy: nicolas.hardy.hernandez@gmail.com
- Dimitris Korobilis: dikorobilis@googlemail.com
- Or open an issue on this GitHub repository
