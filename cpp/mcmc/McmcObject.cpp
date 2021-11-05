#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "McmcObject.h"
#include "Household.h"

using namespace std;


//----------------------------------------------------------------------
// Distributions
//----------------------------------------------------------------------
// Probability density distributions
double logdlnorm(double x, double mean, double sd) {
    return -log(x) - log(sd) - 0.918938533 - pow(log(x)-mean, 2) / (2.0 * pow(sd, 2));
}

double logdexp(double x, double rate) {
    return log(rate) - rate * x;
}

// Random generators
double rnorm(std::mt19937_64& gen, double sd) {
    std::normal_distribution<double> dist(0.0, sd);
    return dist(gen);
}

/*double runif(std::mt19937_64& gen) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}*/

//double rlnorm(std::mt19937_64& gen, double mean, double var) {
//    std::lognormal_distribution<double> dist(mean, var);
//    return dist(gen);
//}

//----------------------------------------------------------------------
// MCMC class
//----------------------------------------------------------------------
// Constructor
McmcObject::McmcObject() {
    m_gen.seed(123);
    m_globalLogLik = 0.0;
    m_iterations = 0;
    m_iterTimeInfection = 0;
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;
    m_rateForRandomWalk.resize(0);
    m_numberOfMoveAccepted.resize(0);
    m_numberOfMoveProposed.resize(0);
    m_parameter.resize(0);
    m_selectedParam.resize(0);
    m_hhLogLik.resize(0);
    std::vector<Household> m_data;
    m_data.resize(0);
    m_nHH = 0;
    m_mainHHSize = 2.0;
    m_sdLNormInfPrior = 1.0;
    m_sdLNormSPrior = 1.0;
    m_maxPCRDetectability = 10.0;
}

McmcObject::McmcObject(
                       unsigned seed,
                       int nIterations,
                       std::vector<Household> data,
                       std::vector<double> parameter,
                       std::vector<int> selectedParameter,
                       std::vector<double> rateForRandomWalk,
                       int nIterTimeInfection, 
                       double mainHHSize, 
                       double sdLNormInfPrior,
                       double sdLNormSPrior,
                       double maxPCRDetectability
                       ) {
    // Default values
    m_globalLogLik= 0.0;
    m_numberOfMoveAccepted.resize(0);
    m_numberOfMoveProposed.resize(0);
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;

    // Copy input vectors
    m_gen.seed(seed);
    m_hhLogLik.resize(data.size());
    m_iterations = nIterations;
    m_iterTimeInfection = nIterTimeInfection;
    m_mainHHSize = mainHHSize;
    m_rateForRandomWalk = rateForRandomWalk;
    m_sdLNormInfPrior = sdLNormInfPrior;
    m_sdLNormSPrior = sdLNormSPrior;
    m_maxPCRDetectability = maxPCRDetectability;
    m_parameter = parameter;
    m_selectedParam = selectedParameter;
    m_data = data;
    m_nHH = data.size();
}

// Initialize parameter values
void McmcObject::initial_param_values() {
    
    for (size_t parID=0; parID < m_parameter.size(); parID++) {
        if ( m_selectedParam[parID] ) {
            
            switch (parID) {
            case 0: // alpha - Uniform(0,1)
                m_parameter[parID] = runif(m_gen); 
                break;

            case 1: // beta - Uniform(0,5)
                m_parameter[parID] = runif(m_gen, 0.0, 5.0); 
                break;

            case 2: // delta - Uniform(-3,3)
                m_parameter[parID] = runif(m_gen, -3.0, 3.0);
                break;

            case 3: // rSAV - LogNormal(0,sd)
                m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormSPrior); 
                break;

            case 4: // rSAVI - LogNormal(0,sd)
                m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormSPrior); 
                break;

            case 5: // rSAI - LogNormal(0,sd)
                m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormSPrior);
                break;

            case 6: // rSC - LogNormal(0,sd)
                m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormSPrior);
                break;

            case 7: // rSCI - LogNormal(0,sd)
                m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormSPrior);
                break;

            case 8: // rInfVac - LogNormal(0,sd)
                m_parameter[parID] = rlnorm(m_gen, 0.0, m_sdLNormInfPrior);
                break;

            case 9: // rAsymptomatic - LogNormal(0,1.0)
                m_parameter[parID]  = rlnorm(m_gen, 0.0, 1.0);
                break;
            }
        }
    }
    
}

// Modification of attributes
void McmcObject::resetMoves() {
    m_numberOfMoveProposed.clear();
    m_numberOfMoveProposed.resize(m_parameter.size());
    m_numberOfMoveAccepted.clear();
    m_numberOfMoveAccepted.resize(m_parameter.size());
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;
}

// Initialize infection time
void McmcObject::initialize_inf_time() {
    
    int ind;
    size_t house;
    for (house=0; house<m_nHH; house++) {

        // Initialize infection time
        int nInd = m_data[house].getSize();
        for (ind=0; ind<nInd; ind++) {
            m_data[house].initialInfTime(ind, m_maxPCRDetectability, m_gen);
        }

        // Compute person-to-person transmission rate within household
        int display = 0;
        m_data[house].compute_lambdas(display);

    }
}

// Initialize infection time
void McmcObject::initial_log_lik() {
    
    size_t house;
    for (house=0; house < m_nHH; house++) {
        m_hhLogLik[house] = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, m_mainHHSize);
        m_globalLogLik += m_hhLogLik[house];
    }
}


// Update parameter value in the mcmc chain
void McmcObject::update_parameter(int parID) {

    // Random walk
	double oldValue = m_parameter[parID];
	double newValue(0.0);
    double logRatioProposal = 0.0;
	if (parID != 2) { // Lognormal random walk for all positive parameters 
		newValue = oldValue * exp(m_rateForRandomWalk[parID] * rnorm(m_gen, 1.0));
        logRatioProposal = log(newValue) - log(oldValue);
	} else { // Normal random walk for delta parameter 
		newValue = oldValue + rnorm(m_gen, m_rateForRandomWalk[parID]);
	}


    // Log ratio of priors 
    double logRatioPrior = 0.0;
    switch (parID) {
    	case 0: // alpha - Uniform(0,1)
    		if (newValue > 1) logRatioPrior = log(0); 
    		break;

        case 1: // beta - Uniform(0,1)
        	if (newValue > 5) logRatioPrior = log(0); 
        	break;

        case 2: // delta - Uniform(-3,3)
        	if (newValue > 3 || newValue < -3) logRatioPrior = log(0);
        	break;

        case 3: // rSAV - LogNormal(0,sd)
        	logRatioPrior = logdlnorm(newValue, 0.0, m_sdLNormSPrior) - logdlnorm(oldValue, 0.0, m_sdLNormSPrior);
        	break;

        case 4: // rSAVI - LogNormal(0,sd)
            logRatioPrior = logdlnorm(newValue, 0.0, m_sdLNormSPrior) - logdlnorm(oldValue, 0.0, m_sdLNormSPrior);
            break;

        case 5: // rSAI - LogNormal(0,sd)
            logRatioPrior = logdlnorm(newValue, 0.0, m_sdLNormSPrior) - logdlnorm(oldValue, 0.0, m_sdLNormSPrior);
            break;

        case 6: // rSC - LogNormal(0,sd)
            logRatioPrior = logdlnorm(newValue, 0.0, m_sdLNormSPrior) - logdlnorm(oldValue, 0.0, m_sdLNormSPrior);
            break;

        case 7: // rSCI - LogNormal(0,sd)
            logRatioPrior = logdlnorm(newValue, 0.0, m_sdLNormSPrior) - logdlnorm(oldValue, 0.0, m_sdLNormSPrior);
            break;

        case 8: // rInfVac - LogNormal(0,sd)
        	logRatioPrior = logdlnorm(newValue, 0.0, m_sdLNormInfPrior) - logdlnorm(oldValue, 0.0, m_sdLNormInfPrior);
            break;

        case 9: // rAsymptomatic - LogNormal(0,1.0)
        	logRatioPrior = logdlnorm(newValue, 0.0, 1.0) - logdlnorm(oldValue, 0.0, 1.0);
            break;
    }


	// Compute new log likelihood
    m_parameter[parID] = newValue;
    std::vector<std::string> parameterNames = {"alpha", "beta", "delta", "rSAV", "rSAVI", "rSAI", "rSC", "rSCI", "rInfVac", "rAsymptomatic"};
    std::vector<double> newLogLikHH(m_nHH, 0.0);
    double newLogLikGlobal(0.0);
    size_t house;
    for (house=0; house < m_nHH; house++) {
        int display = 0;
        newLogLikHH[house] = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, m_mainHHSize, display);        
        newLogLikGlobal += newLogLikHH[house];
    }

	// Update log likelihood
	double Q = newLogLikGlobal - m_globalLogLik + logRatioProposal + logRatioPrior;

	if ( log(runif(m_gen)) < Q )
	{
		m_globalLogLik = newLogLikGlobal;
		m_hhLogLik = newLogLikHH;
        m_numberOfMoveAccepted[parID]++;
	}
	else {
		m_parameter[parID] = oldValue;
	}
	
    m_numberOfMoveProposed[parID]++;
}



// Update augmented data
void McmcObject::update_augmented_data() {

    // Loop over houses
    size_t house, ind;
    for (house = 0; house<m_data.size(); ++house) {
        std::vector<int> infected = m_data[house].getInfectedIndex();

        // Loop over infected individuals
        for (ind = 0; ind < infected.size(); ind++) {

            int inf = infected[ind];
            int display = 0;

            // Independent sampler from the model distribution
            double oldValue = m_data[house].getSpInfTime(inf);
            double newValue = m_data[house].newInfTime(inf, m_maxPCRDetectability, m_gen);

            // Generation ratio
            double oldProposal = m_data[house].pIncub(inf, oldValue, m_maxPCRDetectability);
            double newProposal = m_data[house].pIncub(inf, newValue, m_maxPCRDetectability);
            double logRatioProposal = log(oldProposal)-log(newProposal);

            // Update infection time and person-to-person transmission rates
            m_data[house].setInfTime(inf, newValue);
            m_data[house].update_lambdas(inf, display);   
            
            // Log likelihood
            double newLogLikOfHH = m_data[house].compute_log_lik(m_parameter, m_selectedParam, m_maxPCRDetectability, m_mainHHSize, display);
            double currentLogLik = m_hhLogLik[house];
            double differenceLogLik = newLogLikOfHH - currentLogLik;

            if( log(runif(m_gen)) < differenceLogLik+logRatioProposal ) {
                m_globalLogLik += differenceLogLik;
                m_hhLogLik[house] = newLogLikOfHH;
                m_acceptedMoveData++;

            } else {                
                m_data[house].setInfTime(inf, oldValue);
                m_data[house].update_lambdas(inf);
                
            }

            m_proposedMoveData++;
        }
    }
}

