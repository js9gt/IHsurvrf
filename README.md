# IHsurvrf

IHsurvrf is an R package developed by Jane She for her doctoral dissertation work under the guidance of Michael R Kosorok at the Department of Biostatistics at UNC Chapel Hill.

## Description

The package is designed for estimating optimal dynamic treatment regimes for survival outcomes in multi-stage data. More specifically, it was designed for the purposes of indefinite horizon settings.

The method works through maximizing the truncated mean criterion, and operates through fitting random forest models, with the optional use of multiple strata to divide observations.

For more details, please refer to the original paper.

## Installation:

Method 1. The package can be installed in R directly through github using the following:

``` r
devtools::install_github("js9gt/IHsurvrf/Scripts/IHsurvrf")
```

Then, Load the package in R by running: `library(IHsurvrf)`

Method 2: The package can also be installed after downloading the tar.gz file, located [here](github.com/js9gt/IHsurvrf/Scripts/IHsurvrf_0.1.0.tar.gz). The following steps detail how to install the tar.gz file:

1.  Load R version 4.4.0 in Linux using the bash command: `module load r/4.4.0`
2.  Install the tar.gz file in R by running: `install.packages('[location of .tar.gz file]/IHsurvrf_0.1.0.tar.gz', repos = NULL, type = 'source')`
3.  Load the package in R by running: `library(IHsurvrf)`
4.  The following packages should be automatically loaded once step 3 is run: `survival` `dplyr` `tidyr`

## Data Requirements:

The stage labels separating decision points must use "\_", as does the input stageLabel argument.

When using more than 1 strata, the package also requires a "cumulative.time"" variable denoting how much time a patient has experienced in the study.

Additionally, when using more than 1 strata, the censoring indicator must be defined in the dataset as "delta". Otherwise when just using a single strata, it can have any name.

Finally, the package requires a variable representing patient IDs to be labeled as "subj.id".

## Simulation Studies:

R code to generate multistage data and evaluate the use of our method can be found [here](github.com/js9gt/IHsurvrf/tree/main/Scripts/Data%20Simulations).

## Data Application Code:

It should first be noted, that to our knowledge, there are no current value function estimators for the survival setting when the number of stages is large. This is an open research area. Therefore, we demonstrate the use of our method using standard IPCW with 3 stages maximum. Though this does not display the strength of our method, which was intended for a large number of stages, it demonstrates a proof-of-concept of our method. The example code can be found [here](https://github.com/js9gt/RL_IBS).
