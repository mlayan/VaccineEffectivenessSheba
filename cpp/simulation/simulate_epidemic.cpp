#include <Rcpp.h>
#include "simulate_epidemic_fun.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame hhEpidemic(
	DataFrame H,
	double beta,
	double alpha,
	double delta, 
	double rSAVI,
	double rSAV, 
	double rSAI,
	double rSC,
	double rSCI,
	double rInfVac,
	NumericVector pAsymptomatic,
	double mainHHSize, 
	double dt
  )
  {
  // Simulate epidemic within a household

  // Data info
  int hhsize = H.nrows();
  IntegerVector isolated = H["isolation"];
  IntegerVector vaccinated = H["vaccinated"];
  IntegerVector adult = H["adult"];
  NumericVector studyPeriods = H["studyPeriod"];

  // Initialize transmission vectors
  IntegerVector infectionStatus = H["infectionStatus"];
  NumericVector dds_temp = H["dds"]; // Symptom onset for symptomatic cases or detection date for asymptomatic cases
  NumericVector dds = clone(dds_temp); 
  NumericVector ddi = rep(1000.0, hhsize); // Infection date

  // Susceptible and index cases at the start of the epidemic
  IntegerVector infectors(0);
  IntegerVector sus(0);
  double lastDate = max(studyPeriods); // End of follow-up
  double firstInfection = 1000.0;

  for (int ind=0; ind<hhsize; ind++) {

  	if (dds[ind] != 1000.0) { // Index cases
		infectors.push_back(ind);
    		if (infectionStatus[ind] == 1) { // Symptomatic index case
    			ddi[ind] = rIncub(dds[ind]);
    		} else { 			// Asymptomatic index case
    			ddi[ind] = rInfectionAsymptomatic(dds[ind]);
    		}

    		// First infection time
    		firstInfection = min(firstInfection, ddi[ind]);
    	
    	} else { // At risk members
    		sus.push_back(ind);
    	}
  }

  // Agent based epidemic within household
  for (int n = 0; n <= (lastDate-firstInfection)/dt; ++n) { 
    
    // Current time
    double curr_time = firstInfection + dt * n;  
    
    // Draw new infections
    IntegerVector newInfected;
    NumericVector infected = runif(sus.size());
    
    // Force of infection of at risk individuals at curr_time 
    for (int s = 0; s < sus.size(); ++s) {
      
      int ind = sus[s]; // Individual

      // Probability of getting infected
      NumericVector FOIS = foi(
        curr_time, 
        dt, 
        studyPeriods[ind],
        dds[infectors], 
        ddi[infectors], 
        vaccinated[infectors], 
        infectionStatus[infectors], 
        beta, 
        alpha, 
        delta, 
        rInfVac, 
        hhsize,
        mainHHSize
        );
      
      // Update foi according to vaccination, isolation and age 
      NumericVector fois = clone(FOIS);

  		if (fois.size() > 1) { // fois[0] corresponds to the infection within the community
  			if ( vaccinated[ind] == 1 && isolated[ind] > 0 && adult[ind] == 1 ) {
			// Relative susceptibility of isolated vaccinated adults/teenagers 
  				for (int i = 1; i<fois.size(); ++i) fois[i] *= rSAVI;

  			} else if ( vaccinated[ind] == 1 && isolated[ind] < 1 && adult[ind] == 1) { 
			// Relative susceptibility of unisolated vaccinated adults/teenagers
  				for (int i = 1; i<fois.size(); ++i) fois[i] *= rSAV;

  			} else if ( vaccinated[ind] == 0 && isolated[ind] > 0 && adult[ind] == 1 ) { 
			// Relative susceptibility of isolated unvaccinated adults/teenagers
  				for (int i = 1; i<fois.size(); ++i) fois[i] *= rSAI;

  			} else if ( vaccinated[ind] == 0 && isolated[ind] > 0 && adult[ind] == 0 ) { 
			// Relative susceptibility of isolated children under 12 y.o.
  				for (int i = 1; i<fois.size(); ++i) fois[i] *= rSCI;

  			} else if ( vaccinated[ind] == 0 && isolated[ind] < 1 && adult[ind] == 0) { 
			// Relative susceptibility of unisolated children under 12 y.o.
  				for (int i = 1; i<fois.size(); ++i) fois[i] *= rSC;

  			}
  		}

      
      double pInfection = 1-exp( -sum(fois) ); // Probability of infection at curr_time
      double asymptomatic = runif(1)[0]; 

      if ( infected[s] < pInfection) { // Ind is infected
        // Update ddi and dds
        ddi[ind] = curr_time + dt*runif(1)[0];

        if (asymptomatic < pAsymptomatic[hhsize-2]) { // Asymptomatic case
        	infectionStatus[ind] = 2;
        	dds[ind] = ddi[ind] + detectionPeriod();
        } else {					// Symptomatic case
        	infectionStatus[ind] = 1;
        	dds[ind] = ddi[ind] + incubPeriod();
        }
             
        // Update new infected individuals
        newInfected.push_back(ind);
      }
    }
    
    
    // Update sus and infectors
    for (int index = 0; index < newInfected.size(); ++index) {
    	infectors.push_back(newInfected[index]);
    	for (int i = 0; i < sus.size(); ++i) {
    		if (newInfected[index] == sus[i]) sus.erase(i);
        }
    }

  }
  
  // Output in R data frame format
  return DataFrame::create(
  	_("indid") = H["indid"],
  	_("hhid") = H["hhid"],
  	_("index") = H["index"],
  	_("hhsize") = H["hhsize"],
  	_("ddi") = ddi,
  	_("dds") = dds,
  	_("infectionStatus") = infectionStatus,
  	_("vaccinated") = vaccinated,
	_("isolation") = isolated,
  	_("studyPeriod") = studyPeriods,
  	_("adult") = adult
  	);
}




