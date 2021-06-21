#ifndef SIMULATE_EPIDEMIC_FUN_PKG__H
#define SIMULATE_EPIDEMIC_FUN_PKG__H

#include <Rcpp.h>
using namespace Rcpp;

extern double mIncub;
extern double sdIncub;
extern double maxPCRDetectability;

RcppExport double incubPeriod();
RcppExport double detectionPeriod();
RcppExport double rIncub(double d);
RcppExport double rInfectionAsymptomatic(double d);
RcppExport NumericVector foi(
	double t, 
	double dt,
	double lastDate,  
	NumericVector symptomOnset, 
	NumericVector infection, 
	IntegerVector vaccinatedInfectors,
	IntegerVector infStatus, 
	double beta, 
	double alpha, 
	double delta,
	double rInfVac, 
	double hhsize,
	double mainHHSize
    );

#endif