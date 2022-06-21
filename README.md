# Impact of BNT162b2 vaccination and isolation on SARS-CoV-2 transmission in Israeli households: an observational study

The full article is available [here](https://doi.org/10.1093/aje/kwac042).

## Dependencies

* A C++ compiler, like gcc
* The boost C++ libraries
* R

The code has been tested using g++ (8.4.0, 9.2.0), Boost (1.72.0), and R (4.0.2) on Linux CentOS 6.10.

## Compilation

Set the variable CC in the Makefile to point to the compiler installed on your system

In the `cpp/mcmc/` directory, compile the code by typing (in your terminal): `make`

The compilation takes a few seconds and creates the executable `program.out`.

## Usage

### Data

The space-separated file `data/2021_05_14_full_database_2doses.txt` contains the original data used in the article. It has 10 columns (no header, one row per individual) indicating:

* the individual id 
* the household id 
* the household size
* the date of symptom onset for symptomatic cases or the date of detection for asymptomatic cases (1000: not infected) 
* the infection status (0: not infected; 1: symptomatic; 2: asymptomatic; 3: symptomatic with missing symptom onset)
* the vaccination status according to the definition of effective vaccination (-1: infected in the preceding year; 0: no; 1: yes) 
* the end of the follow-up
* the age of the individual (0: child under 12 y.o.; 1: adult/teenager above 12 y.o.)
* the index cases (0: index case; 1: household contact)
* the isolation status (0: no isolation; 1: partial isolation; 2: complete isolation)

Additional data files were used in the sensitivity analysis (Figure 4 and Supplementary Materials):

* `2021_05_14_full_database_1dose.txt` contains the same households than the principal anlaysis and vaccination is assumed effective 15 days after the 1st dose 
* `2021_05_14_1PCR_2doses.txt` contains the households where all negative contacts performed at least one PCR test in the 10 days following the detection of the index case
* `2021_05_14_2PCR_2doses.txt` contains the households where all negative contacts performed at least two PCR tests in the 10 days following the detection of the index case
* `2021_05_14_no_early_vaccination_2doses.txt` does not contains the households where the index case was vaccinated but got infected before the vaccine was considered effective (>7 days after the 2nd dose)

### Launching the analysis

Use the `launch_mcmc.R` script in the `cpp/mmc/` folder to set up and launch the analyses.

Running the script as provided takes about 20 minutes dependending on the number of participants that is analyzed. 
The output is written to a space-delimited file in the `results/` folder. Its name contains the arguments passed to the program in the following order: database, vaccination definition, log-sd of the relative infectivity/relative susceptibility prior, relative infectivity of asymptomatic cases and chain id. Each file contains one MCMC chain. 

### Simulating household epidemics

The script `R/3.simulation.R` allows one to simulate 2,000 household epidemics based on the original data `data/2021_05_14_full_database_2doses.txt` and the Rcpp scripts in the `cpp/simulation` directory. Running the script as provided takes about 40 minutes. 

This script depends on two R libraries: 1) tidyverse and 2) Rcpp.

The script `R/4.model_fit.R` allows one to plot the fit of the model that is stored in the `figures/` directory.

This script depends on two R libraries: 1) tidyverse and 2) binom.

### Visualizing the results

The scripts `R/1.output_analysis.R` and `R/2.sensitivity_analysis.R` allows one to visualize parameter estimates in the baseline scenario and in the sensitivity analysis, respectively. All figures and tables are generated in the `figures/` and `tables/` directories. 

These scripts depend on two R libraries: 1) tidyverse and 2) gridExtra. 

## Erratum
There was an error in the likelihood of the inference model at line 575 of `Household.cpp` where the probability of infection was written in a discrete way, $1-exp⁡(-λ(t_{infection}))$, 
while our model is in continuous time and thus the probability of infection is equal to $λ(t_{infection})$. The code has been corrected and all output files available in the repository have been updated. For more details on the impact of this error on the estimates, we refer to the published erratum.

