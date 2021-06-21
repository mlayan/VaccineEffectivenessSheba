# Impact of BNT162b2 vaccination and isolation on SARS-CoV-2 transmission in Israeli households: an observational study

The full article is available [here]().

## Dependencies

* A C++ compiler, like gcc
* The boost C++ libraries
* R

The code has been tested using g++ (7.2.0), Boost (1.72.0), and R (4.0.2) on Linux CentOS 6.10.

## Compilation

Set the variable CC in the Makefile to point to the compiler installed on your system

Compile the code by typing (in your terminal): `make`

The compilation takes a few seconds and creates the executable program.out.

## Usage

### Data

The space-separated file `data/2021_05_14_full_database_2doses.txt` contains the original data used for the article. It has 10 columns (no header, one row per individual) indicating:

* the individual id 
* the household id 
* the household size
* the date of symptom onset or the date of detection for asymptomatic cases 
* the infection status (0: not infected, 1: symptomatic, 2: asymptomatic, 3: symptomatic with missing symptom onset)
* the vaccination status according to the definition of effective vaccination (0: no, 1: yes) 
* the end of the follow-up
* the age of the individual (0: child, 1: adult > 12 y.o.)
* the index cases (0: index case, 1: household contact)
* the isolation status (0: no isolation, 1: partial isolation, 2: complete isolation)


Additional data files were used in the sensitivity analysis (Figure 4 and Supplementary Materials):

* `2021_05_14_full_database_1dose.txt` contains the same households than the principal anlaysis and vaccination is assumed effective >15 days after the 1st dose 
* `2021_05_14_known_outcome_2doses.txt` contains the households where all contacts performed a PCR test in the 10 days following the detection of the index case
* `2021_05_14_full_database_2doses_strict.txt` does not contains the households where the index case was vaccinated but got infected before the vaccine was considered effective (>7 days after the 2nd dose)

### Launching the analysis

Run the R script launchLocal.R using your preferred R environment. This script allows one to set all the different options to run the analyses (e.g. the number of MCMC steps or the choice of change points for R0).
Running the script as provided takes about 20 minutes.
The output will be written to there files `results/mcmc_1.txt`, `results/mcmc_2.txt` and `results/mcmc_3.txt` that contain one MCMC chain each.
Note that the analysis can also be performed without using R, by simply passing all the appropriate options as input arguments to covid.exe.

### Simulating household epidemics


### Visualizing the results

The scripts 1.output_analysis.R, 3.model_adequacy.R and 4..R allows one to visualize the estimations. 
This scripts depends on two R libraries: 1) tidyverse; and 2) gridExtra. 
