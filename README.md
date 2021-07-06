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

The compilation takes a few seconds and creates the executable `program.out`.

## Usage

### Data

The space-separated file `data/2021_05_14_full_database_2doses.txt` contains the original data used in the article. It has 10 columns (no header, one row per individual) indicating:

* the individual id 
* the household id 
* the household size
* the date of symptom onset for symptomatic cases or the date of detection for asymptomatic cases (1000: not infected) 
* the infection status (0: not infected, 1: symptomatic, 2: asymptomatic, 3: symptomatic with missing symptom onset)
* the vaccination status according to the definition of effective vaccination (0: no, 1: yes) 
* the end of the follow-up
* the age of the individual (0: child under 12 y.o., 1: adult/teenager above 12 y.o.)
* the index cases (0: index case, 1: household contact)
* the isolation status (0: no isolation, 1: partial isolation, 2: complete isolation)


Additional data files were used in the sensitivity analysis (Figure 4 and Supplementary Materials):

* `2021_05_14_full_database_1dose.txt` contains the same households than the principal anlaysis and vaccination is assumed effective 15 days after the 1st dose 
* `2021_05_14_known_outcome_2doses.txt` contains the households where all contacts performed a PCR test in the 10 days following the detection of the index case
* `2021_05_14_full_database_2doses_strict.txt` does not contains the households where the index case was vaccinated but got infected before the vaccine was considered effective (>7 days after the 2nd dose)

### Launching the analysis

Pass the following options as input arguments to `program.out` in your terminal:

* Length of the MCMC chain
* ID of the MCMC chain (1; 2; 3)
* log-sd of the log-normal prior distribution of the relative infectivity parameter
* log-sd of the log-normal prior distribution of the relative susceptibility parameters
* Relative infectivity of asymptomatic cases compared to symptomatic cases
* Effective vaccination definition (1dose; 2doses)
* Database (full_database; know_outcome; strict)

Running the script as provided takes about 20 minutes.
The output will be written to a space-delimited file in the `results/` folder. Its name contains the arguments passed to the program in the following order: database,  vaccination definition, log-sd of the relative infectivity prior, log-sd of the relative susceptibility prior, relative infectivity of asymptomatic cases and chain id. Each file contains one MCMC chain. 

### Simulating household epidemics


### Visualizing the results

The scripts `1.output_analysis.R`, `3.model_adequacy.R` and `4.sensitivity_analysis.R` allows one to visualize the estimations. 
These scripts depend on two R libraries: 1) tidyverse and 2) gridExtra. 
