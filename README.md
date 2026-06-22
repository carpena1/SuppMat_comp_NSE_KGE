# Supplementary Materials

## Theoretical connections between the Nash–Sutcliffe efficiency and the Kling–Gupta Efficiency for environmental model evaluation

By Axel Ritter<sup>1</sup> and Rafael Muñoz-Carpena<sup>2,\*</sup>

*Under review in Journal of Hydrology, June 2026*

This repository contains the supplementary files with the statistics and metrics for all the datasets used for model evaluation in the manuscript (Table 1).

### Datasets 

**Datasets description**

The summary model efficiency sttistics are summarised in Table 1 from the manuscript shown below.

**Table 1. Model testing datasets used.**

| Dataset Nr./ Dataset Ref. | Predicted variables <sup>a</sup> | N <sup>b</sup> | NSE range <sup>c</sup> [KSE range] | Reference |
|---|---|---|---|---|
| **Testing** | | | | |
| 1/FP-WQ | Ground water [F<sup>-</sup>], [Cl<sup>-</sup>], [N-NO<sub>3</sub><sup>-</sup>], [P-PO<sub>4</sub><sup>3-</sup>], [TP] (-) | 50 | 0.215 – 0.994 [-10<sup>14</sup> – 0.836] | Muñoz-Carpena et al. (2005) |
| 2/LOX-EC | Ground water electrical conductivity (-) | 7 | 0.569 – 1.000 [-0.532 – 1.000] | Kaplan and Muñoz-Carpena (2014) |
| 3/LOX-SM | Soil moisture (-) | 12 | 0.629 – 0.983 [-10<sup>2</sup> – 0.931] | Kaplan and Muñoz-Carpena (2011) |
| 4/LOX-WT | Water table elevation (-) | 12 | 0.778 – 1.000 [-1.794 – 1.000] | Kaplan et al. (2010) |
| 5/OKA-NDVI | Normalized Difference Vegetation Index, NDVI (-) | 240 | 0.508 – 0.931 [-0.261 – 0.623] | Campo-Bescós et al. (2013) |
| 6/SH-WQ | Ground water [N-NO<sub>3</sub><sup>-</sup>] (mg/l) | 18 | 0.462 – 0.985 [0.529 – 0.968] | Ritter et al. (2007) |
| 7/WDPT– 2T | Water droplet penetration time (s) | 80 | 0.546 – 0.967 [0.631 – 0.976] | Regalado and Ritter (2009a) |
| 8/WDPT– 3T | Water droplet penetration time (s) | 80 | 0.549 – 0.991 [0.633 – 0.994] | Regalado and Ritter (2009a) |
| 9/MED | Contact angle (º) | 40 | 0.893 – 0.994 [0.923 – 0.996] | Regalado and Ritter (2009b) |
| **Verification** | | | | |
| 10/CAMELS | Runoff (mm d<sup>-1</sup>) | 613 | 0.207 – 0.922 [-0.022 – 0.933] | Newman et al. (2015) |

<sup>a</sup> (-) indicates dimensionless standardized data series, which implies mean of observation close to zero, and thus values with alternating sign in the dataset; <sup>b</sup> N is the number of model evaluation in each study; <sup>c</sup> Performance metrics were computed using the FITEVAL software (Ritter and Muñoz-Carpena (2013)) option that allows for discarding repeated paired values.

**Download Table:**
The extended Table containing the metrics for each set in the dataset used in the study figures and Tables can be downloaded here **[Datasets_1-10_statistics.csv](Datasets_1-10_statistics.csv)**

### Analysis tools
FITEVAL (https://abe.ufl.edu/carpena/software/fiteval/) Matlab scripts were used for analysis of the individual {obs_i, pred_1} results in each study of the datasets in Table 1.
