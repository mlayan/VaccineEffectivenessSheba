#ifndef MCMCOBJECT__H
#define MCMCOBJECT__H

#include <random>
#include <vector>
#include "Household.h"

// Useful distributions
double logdlnorm(double x, double mean, double var);
double runif(std::mt19937_64& gen, double lower_value, double upper_value);
double rnorm(std::mt19937_64& gen, double mean, double var);
double rlnorm(std::mt19937_64& gen, double mean, double var);

// Class
class McmcObject
{
public:
	McmcObject();
	McmcObject(
            size_t seed,
            int nIterations,
            std::vector<Household> data,
            std::vector<double> parameter,
            std::vector<int> selectedParameter,
            std::vector<double> rateForRandomWalk,
            double mainHHsize = 2.0,
            double sdLNormInfPrior = 1.0,
            double sdLNormSPrior = 1.0
            );
	~McmcObject() {};

	int getNumbHH() const { return m_data.size(); }
	int iteration() const { return m_iterations; };
	int proposedMove(int index) const {return m_numberOfMoveProposed[index]; };
	int acceptedMove(int index) const {return m_numberOfMoveAccepted[index]; };
	int proposedMoveData() const {return m_proposedMoveData; };
	int acceptedMoveData() const {return m_acceptedMoveData; };
	double globalLogLik() const { return m_globalLogLik; };
	double hhLogLik(int index) const { return m_hhLogLik[index]; };
	std::vector<double> hhLogLik() const { return m_hhLogLik; };
	double rateRandomWalk(int index) const { return m_rateForRandomWalk[index]; };
	double parameter(int index) const { return m_parameter[index]; };
	size_t nParameters() const { return m_parameter.size(); };

	void resetMoves();
	void initialize_inf_time();
	void initial_log_lik();
	void update_parameter(int parID);
	void update_augmented_data();

private:
	int m_iterations;
	std::mt19937_64 m_gen;
	double m_globalLogLik;
	std::vector<double> m_hhLogLik;
	int m_acceptedMoveData;
	int m_proposedMoveData;
	std::vector<int> m_numberOfMoveAccepted;
	std::vector<int> m_numberOfMoveProposed;
	std::vector<double> m_parameter;
	std::vector<int>  m_selectedParam;
	std::vector<double> m_rateForRandomWalk;
	double m_sdLNormInfPrior;
	double m_sdLNormSPrior;
	std::vector<Household> m_data;
	size_t m_nHH;
	double m_mainHHSize;
};

#endif
