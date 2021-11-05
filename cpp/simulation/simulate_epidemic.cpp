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
	double rSAV,
  double rSAVI,
  double rSAI,
  double rSC,
  double rSCI,
	double rInfVac,
  double rAsymptomatic,
	double pAsymptomatic,
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
  double lastDate = max(studyPeriods);

  // Initialize transmission vectors
  IntegerVector infectionStatus = H["infectionStatus"];
  NumericVector dds_temp = H["dds"];
  NumericVector dds = clone(dds_temp);
  NumericVector ddi = rep(1000.0, hhsize);

  // Susceptible and index cases at the start of the epidemic
  IntegerVector infectors(0);
  IntegerVector sus(0);
  double firstInfection = 1000.0;

  for (int ind=0; ind<hhsize; ind++) {

    if (dds[ind] != 1000.0) { // Index cases

      infectors.push_back(ind);
      if (infectionStatus[ind] == 1) { // Symptomatic cases
        ddi[ind] = rIncub(dds[ind]);
      } else { // Asymptomatic cases
        ddi[ind] = rInfectionAsymptomatic(dds[ind]);
      }

      // First infection time
      // Corresponds to the start of the in silico follow-up of the household
      firstInfection = min(firstInfection, ddi[ind]); 

    } else { // Susceptibles
      sus.push_back(ind);
    }
  }

  // Agent based epidemic within household
  for (int n = 0; n <= (lastDate-firstInfection)/dt; ++n) {
    
    // Current time
    double curr_time = ddi[0] + dt * n;

    // Draw new infections
    IntegerVector newInfected;
    NumericVector infected = runif(sus.size());  

    // Vector of fois in the community and partial foi from the 
    // infected households members at curr_time
    NumericVector FOIS = foi(
      curr_time, 
      dt, 
      lastDate,
      dds[infectors], 
      ddi[infectors], 
      adult[infectors],
      vaccinated[infectors], 
      infectionStatus[infectors], 
      beta, 
      alpha, 
      delta, 
      rInfVac,
      rAsymptomatic,
      hhsize,
      mainHHSize
      );

    // Loop on susceptibles 
    for (int s = 0; s < sus.size(); ++s) {
      // Individual
      int ind = sus[s];
      
      // Update foi according to age
      NumericVector fois = clone(FOIS);

      if (fois.size() > 1) { // fois[0] corresponds to the infection within the community
        if ( vaccinated[ind] == 1 && isolated[ind] > 0 && adult[ind] == 1 ) { // Relative susceptibility of isolated vaccinated adults 
          for (int i = 1; i<fois.size(); ++i) fois[i] *= rSAVI;

        } else if ( vaccinated[ind] == 1 && isolated[ind] < 1 && adult[ind] == 1) { // Relative susceptibility of unisolated vaccinated adults
          for (int i = 1; i<fois.size(); ++i) fois[i] *= rSAV;

        } else if ( vaccinated[ind] == 0 && isolated[ind] > 0 && adult[ind] == 1 ) { // Relative susceptibility of isolated unvaccinated adults
          for (int i = 1; i<fois.size(); ++i) fois[i] *= rSAI;

        } else if ( vaccinated[ind] == 0 && isolated[ind] > 0 && adult[ind] == 0 ) { // Relative susceptibility of isolated children
          for (int i = 1; i<fois.size(); ++i) fois[i] *= rSCI;

        } else if ( vaccinated[ind] == 0 && isolated[ind] < 1 && adult[ind] == 0) { // Relative susceptibility of unisolated children
          for (int i = 1; i<fois.size(); ++i) fois[i] *= rSC;

        }
      }

      // Is ind infected at curr_time?
      double pInfection = 1-exp( -sum(fois) );
      double asymptomatic = runif(1)[0];

      if (infected[s]<pInfection) {
        // Update ddi and dds
        ddi[ind] = curr_time + dt*runif(1)[0];

        if (asymptomatic < pAsymptomatic) {
        	infectionStatus[ind] = 2;
        	dds[ind] = ddi[ind] + detectionPeriod();
        } else {
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
          if (newInfected[index] == sus[i]) {
            sus.erase(i);
          }
        }
      }
    }
  
  // Output
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




