#include <Rcpp.h>
#include "simulate_epidemic_fun.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


double mIncub = 1.63;
double sdIncub = 0.41;
double maxPCRDetectability = 10.0;


//////////////////////////////////////////////
// [[Rcpp::export]]
double incubPeriod() 
  {
  // Incubation period
  double d = rlnorm(1, mIncub, sdIncub)[0];
  while(d<3 || d >30) d = rlnorm(1, mIncub, sdIncub)[0];
  return d;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double detectionPeriod() {
  return runif(1, 0, maxPCRDetectability)[0];
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double rIncub(double d) 
  {
  // Draw random incubation period for symptomatic cases
  double pIncub = rlnorm(1, mIncub, sdIncub)[0];
  while (pIncub < 3.0 || pIncub > 30.0) {
    pIncub = rlnorm(1, mIncub, sdIncub)[0];  
  }
  
  return d-pIncub;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double rInfectionAsymptomatic(double d) {
	return runif(1, d-maxPCRDetectability, d)[0];
}


//////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector foi(
	double t, 
	double dt, 
	double lastDate,  
	NumericVector symptomOnset, 
	NumericVector infection, 
  IntegerVector age,
	IntegerVector vaccinatedInfectors,
	IntegerVector infStatus, 
	double beta, 
	double alpha, 
	double delta,
  double rInfVac,
  double rAsymptomatic,
	int hhsize, 
  double mainHHSize
  ) 
  {
  // Force of infection from 
  // 	0: community
  // 	1 to #infected: infected household contacts
  NumericVector fois(symptomOnset.size() + 1);
  
  //Infection by the community
  fois[0] = alpha*dt;
  
  // Relative infectivity
  IntegerVector allVaccinated = rep(1, vaccinatedInfectors.size());
  NumericVector relativeInfectivity = ifelse(vaccinatedInfectors == allVaccinated, rInfVac, 1.0);
    
  // Parameters
  double mBeta = 26.1;
  double vBeta = 7;
  double shift = 25.6;
  double shapeBeta = pow(mBeta,2) / vBeta;
  double scaleBeta = vBeta / mBeta;
  double normCons = 1 - R::pgamma(shift - 3, shapeBeta, scaleBeta, true, false);
  
  
  // Infection by infected individuals within the same household
  if (symptomOnset.size() > 0) {
    // NumericVector last = pmin(lastDate, d + 6);
    for (int index = 0; index < symptomOnset.size(); ++index) {
      double k = 0.0;

      if (infStatus[index] == 1) { // Symptomatic infector
      	if (t >= symptomOnset[index] - 3 && symptomOnset[index] < lastDate) {
      		k += ( R::pgamma(shift + (t + dt - symptomOnset[index]), shapeBeta, scaleBeta, true, false) - R::pgamma(shift + (t - symptomOnset[index]), shapeBeta, scaleBeta, true, false) ) / normCons;

        } else if (t >= symptomOnset[index] - 3 && symptomOnset[index] >= lastDate) {
          k += ( R::pgamma(shift + (t + dt - lastDate), shapeBeta, scaleBeta, true, false) - R::pgamma(shift + (t - lastDate), shapeBeta, scaleBeta, true, false) ) / normCons;

        }

      } else { // Asymptomatic infector 
      	if ( t >= (infection[index] + 2) && (infection[index] + 2) < lastDate ) { // Asymptomatic cases are infectious 2 days after their infection
      		k += ( R::pgamma(shift - 3.0 + (t + dt - infection[index] - 2.0), shapeBeta, scaleBeta, true, false) - R::pgamma(shift - 3.0 + (t - infection[index] - 2.0), shapeBeta, scaleBeta, true, false) ) / normCons ;

          k *= rAsymptomatic;
      	}
      }

      fois[index + 1] = (beta * k * relativeInfectivity[index]);
      if (delta != 0.0) fois[index + 1] /= pow(hhsize / mainHHSize, delta);
    } 
  }
  
  return fois;
}
