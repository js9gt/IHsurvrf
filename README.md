IHsurvrf is an R package developed by Jane She for her doctoral dissertation work under the guidance of Michael R Kosorok at the Department of Biostatistics at UNC Chapel Hill.

The package is designed for estimating optimal dynamic treatment regimes for survival outcomes in multi-stage data. More specifically, it was designed for the purposes of indefinite horizon settings.

The method works through maximizing the truncated mean criterion, and operates through fitting random forest models, with the optional use of multiple strata to divide observations.

For more details, please refer to the original paper.

Installation:
The package can be installed directly through github using the following:


Or, the package can be installed after downloading the tar.gz file, located at github.com/js9gt/IHsurvrf/Scripts/IHsurvrf_0.1.0.tar.gz 
1. Load R version 4.4.0 in linux using bash command: module load r/4.4.0
2. install tar.gz file in R by running: install.packages('[location of .tar.gz file]/IHsurvrf_0.1.0.tar.gz. repos = NULL, type='source')
3. load package by in R by running: library(IHsurvrf)

Data Requirements:

Simulation Studies:
R code to generate multistage data and evaluate the use of our method can be found at github/com/js9gt/IHsurvrf/Scripts/'Data Simulations'.

Data Application Code:
It should first be noted, that to our knowledge, there are no current value function estimators for the survival setting when the number of stages is large. This is an open research area.
Therefore, we demonstrate the use of our method using standard IPCW with 3 stages maximum. Though this does not display the strength of our method, which was intended for a large number of stages,
it demonstrates a proof-of-concept of our method. The example code can be found at github.com/js9gt/RL_IBS.
