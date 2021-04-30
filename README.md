# Effects of mRNA vaccine on SARS-CoV-2 transmission in households

The full article is available [here]().

## Dependencies

* A C++ compiler, like gcc
* The boost C++ libraries
* R

The code has been tested using g++ (7.2.0), Boost (1.72.0), and R (4.0.2) on Linux CentOS 6.10.

## Compilation

Set the variable CC in the Makefile to point to the compiler installed on your system

Compile the code by typing (in your terminal):
`make`

The compilation takes a few seconds and creates the executable program.exe.

## Usage

1. Data

The file Data/guyane20200825.txt contains the original data used for the
article. It has four columns (no header, one row per day) indicating:

the number of new hospital admissions
the number of new ICU admissions
the number of occupied general ward beds
the number of occupied ICU beds

The first line corresponds to April 22nd.
The file Data/figures_4c_S5.csv contains the data used to generate Figure 4c
and Supplementary Figure 5.

2. Launching the analysis

Run the R script launchLocal.R using your preferred R environment. This script
allows one to set all the different options to run the analyses (e.g. the number
of MCMC steps or the choice of change points for R0).
Running the script as provided takes about 4 minutes.
The output will be written to the two files Output/mcmc.csv (MCMC chains) and
Output/mcmc_sims.csv (trajectories).
Note that the analysis can also be performed - although not as conveniently -
without using R, by simply passing all the appropriate options as input
arguments to covid.exe.

3. Visualizing the results

The script visualize.R allows one to load the results of the MCMC run and
to visualize the projections. This scripts depends on two R libraries: 1)
tidyverse; and 2) cowplot.
NOTE: the provided files Output/mcmc.csv and Output/mcmc_sims.csv contain
the output of one run of launchLocal.R with the options described in the
script.
