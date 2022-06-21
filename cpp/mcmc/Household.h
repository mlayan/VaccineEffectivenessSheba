#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H
#include <string>
#include <vector>
#include <random>


// Constant parameters
// extern double maxPCRDetectability;
extern double mIncub;
extern double varIncub;
extern double mGamma;
extern double vGamma;
extern double shift;
extern double shape;
extern double scale;

// PDF
double plnorm(double x) ;
double dlnorm(double x, double mean, double var) ;
double pgamma(double x) ;
double dgamma(double x) ;

// Random generators
double runif(std::mt19937_64& gen, double lower_value = 0.0, double upper_value = 1.0) ;
double rlnorm(std::mt19937_64& gen, double mean, double var) ;


//----------------------------------------------------------------------
// Infectivity profile
//----------------------------------------------------------------------
double infectivityProfile(
	double origin,
	double infectionDate,
	int infectionStatus,
	double tinf,
	int studyPeriod
	);
double cumulativeInfectivity(
	double origin,
	double infectionDate,
	int infectionStatus,
	double tinf,
	int studyPeriod
	);


//----------------------------------------------------------------------
// Household class
//----------------------------------------------------------------------
class Household
{
public:
	Household();
	~Household() {};
	// int getStudyPeriod() const { return m_studyPeriod; };
	int getSize() const { return m_size; };

	size_t nInfected() const { return m_confCase.size(); };
	int getSpInfected(int index) const { return m_infected[index] ; };
	std::vector<int> getSpInfected(std::vector<int> index) const;
	std::vector<int> getAllInfected() const { return m_infected; };
	std::vector<int> getInfectedIndex() const { return m_confCase; };

	double getSpInfTime(int index) const { return m_infTime[index]; };
	std::vector<double> getAllInfTime() const { return m_infTime; };

	void addIndividual(int indid, double onsetTime, int isCase, int vaccinationStatus, int age, int studyPeriod, int identifiedIndex, int isolation);
	void setInfTime(int index, double infTime);
	void newHousehold();
	void initialInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen);
	void displayHH(); 
	void compute_lambdas(int display = 0);
	void update_lambdas(int ind, int display = 0);

	double newInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen);
	double pIncub(int index, double infTime, double maxPCRDetectability);

	double compute_log_lik(std::vector<double> parameters, std::vector<int> selectedParam, double maxPCRDetectability, double mainHHSize, int display = 0);
	double log_pInf(int curr, double t0, std::vector<double> parameter, std::vector<int> selectedParam, double mainHHSize, int display = 0);
	double log_S(int curr, double t0, std::vector<double> parameter, std::vector<int> selectedParam, double mainHHSize, int display = 0);

private:
	int m_size;
	int m_notInfected;
	double m_startFollowUp;
	std::vector<int> m_indid;
	std::vector<double> m_onsetTime;
	std::vector<int> m_infected;
	std::vector<int> m_confCase;
	std::vector<int> m_studyPeriod;
	std::vector<int> m_vaccinationStatus;
	std::vector<int> m_age;
	std::vector<int> m_index;
	std::vector<int> m_isolation;
	std::vector<double> m_infTime;
	std::vector<std::vector<double>> m_cumLambda;
	std::vector<std::vector<double>> m_instLambda;
};

#endif