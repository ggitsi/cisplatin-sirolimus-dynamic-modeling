# Cisplatin–Sirolimus Dynamic Modeling

**Computational Prediction of Sirolimus-Mediated Protection Against Cisplatin-Induced Nephrotoxicity**

## Project Overview

This project investigates the potential protective effects of sirolimus against cisplatin-induced nephrotoxicity using a dynamic autophagy–apoptosis model implemented in MATLAB.

The study extends the mathematical model developed by Tavassoly et al. (2015) by incorporating sirolimus-mediated inhibition of the mechanistic target of rapamycin (mTOR). The simulations evaluate how this intervention affects autophagy, cellular stress, and apoptosis in renal proximal tubular cells.

The project also includes sensitivity analysis and computational optimization to investigate the robustness of the predicted protective effects under different model assumptions.

## Mathematical Model

The model is based on a system of ordinary differential equations (ODEs) describing the regulatory interactions between autophagy and apoptosis.

The sirolimus extension introduces an additional inhibition term into the mTOR regulatory input. Its effect is controlled by an exploratory parameter, `kSiro`, which has not been calibrated using experimental data.

Numerical simulations are performed using MATLAB's `ode15s` solver.

## Simulation Workflow

| Step | Description | MATLAB File |
|---|---|---|
| 1 | Baseline single-cell simulation | `Main.m` |
| 2 | Population simulation of 100 virtual cells | `Main_100Cells.m` |
| 3 | Cisplatin dose-response analysis | `Main_CisplatinDoseResponse.m` |
| 4 | Sirolimus-only model validation | `Main_SirolimusOnlyValidation.m` |
| 4B | Sensitivity analysis of the sirolimus inhibition parameter | `Main_kSiroSensitivity.m` |
| 5 | Cisplatin–sirolimus combination sensitivity analysis | `Main_CombinationSensitivity.m` |
| 6 | Mechanistic dynamics of mTOR, autophagy, stress, and apoptosis | `Main_MechanisticDynamics.m` |
| 7 | Nephrotoxicity-constrained treatment optimization | `Main_TreatmentOptimization.m` |
| 8 | Robustness analysis across uncertain inhibition parameters | `Main_kSiroOptimizationRobustness.m` |

The repository also contains `Fun.m`, representing the baseline model equations, and `Fun_Sirolimus.m`, containing the sirolimus-extended model. An alternative Step 6 implementation is provided in `Main_MechanisticDynamics_v2.m`.

## Requirements

- MATLAB with support for the `ode15s` ODE solver.
- Access to the model and simulation scripts included in this repository.

## How to Run

1. Download or clone the repository.
2. Open MATLAB and set the Current Folder to the repository directory.
3. Ensure the required model functions (`Fun.m` and `Fun_Sirolimus.m`) are accessible.
4. Run the corresponding `Main_*.m` functions according to the simulation workflow.
5. Examine the generated figures, CSV result tables, and MATLAB workspace files where applicable.

**Note:** Steps 7 and 8 currently contain a local Windows output-directory path. Update the `outputFolder` variable to a valid directory on your computer before running these scripts.

## Key Simulation Results

The simulations predicted that adding sirolimus to cisplatin treatment could reduce apoptosis in virtual renal proximal tubular cells.

For a cisplatin model input of C = 6 and 10 µM sirolimus, with `kSiro = 0.005`, the simulated apoptosis percentage at 48 hours decreased from **45% to 6%**.

This response was accompanied by increased autophagosome levels, reduced cellular stress, and lower pro-apoptotic BH3 levels.

Sensitivity analysis further showed that the conservative maximum cisplatin stress input satisfying a 20% apoptosis threshold increased from **5.5 to 10.0** at 20 µM sirolimus across the four tested inhibition coefficients.

## Model Limitations

The findings represent computational predictions rather than experimentally validated therapeutic outcomes.

The sirolimus inhibition coefficient has not been calibrated using experimental data. Additionally, the model does not incorporate drug pharmacokinetics, direct sirolimus toxicity, or the potential effects of autophagy modulation on cancer cells.

The cisplatin input C represents a model-based stress input and should not be directly interpreted as a clinical cisplatin dose. The 20% apoptosis threshold is an analytical criterion rather than an established safety limit.

Further experimental validation is necessary before evaluating potential clinical applications.

## Reference

Tavassoly, I., Parmar, J., Shajahan-Haq, A. N., Clarke, R., Baumann, W. T., & Tyson, J. J. (2015). Dynamic modeling of the interaction between autophagy and apoptosis in mammalian cells. *CPT: Pharmacometrics & Systems Pharmacology, 4*(4), 263–272.

https://doi.org/10.1002/psp4.29

## Acknowledgment

This project builds upon the mathematical framework developed by Tavassoly et al. (2015). The sirolimus extension, simulation workflow, sensitivity analysis, and optimization scenarios were implemented as part of an independent computational modeling study.
