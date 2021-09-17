#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include "Household.h"

using namespace std;

// Constant parameters
double maxPCRDetectability = 10;
double mIncub = 1.63;
double sdIncub = 0.5;
double mGamma = 26.1;
double vGamma = 7;
double shift = 25.6;
double shape = pow(mGamma, 2) / vGamma;
double scale = vGamma / mGamma;

//----------------------------------------------------------------------
// Cumulative and probability density functions
//----------------------------------------------------------------------
double plnorm(double x) {
   boost::math::lognormal_distribution<> d(mIncub, sdIncub);
   return cdf(d, x);
}

double dlnorm(double x, double mean, double sd) {
    return exp( - pow(log(x) - mean, 2) / ( 2 * pow(sd,2) ) ) / (x * sd * 2.50662827);
}

double pgamma(double x) {
   boost::math::gamma_distribution<> d(shape, scale);
   return cdf(d, x);
}

double dgamma(double x) {
   boost::math::gamma_distribution<> d(shape, scale);
   return pdf(d, x);
}

// Random generators
double runif(std::mt19937_64& gen, double lower_value, double upper_value) {
    std::uniform_real_distribution<double> dist(lower_value, upper_value);
    return dist(gen);
}

double rlnorm(std::mt19937_64& gen, double mean, double sd) {
    std::lognormal_distribution<double> dist(mean, sd);
    return dist(gen);
}


//----------------------------------------------------------------------
// Infectivity profile
//----------------------------------------------------------------------
double infectivityProfile(
                          double origin,
                          double infectionDate,
                          int infectionStatus,
                          double tinf,
                          int studyPeriod
                          ) {
			  
  // Instantaneous risk of infection
  // The probability density function was estimated by Aschcroft et al., 2020 (DOI:10.4414/smw.2020.20336)
  double out = 0.0;
  double normCons = 1 - pgamma(shift - 3);

  if (infectionStatus == 1) { // Symptomatic infector with a date of symptom onset

    if (origin-3 < studyPeriod && origin-3 < tinf) {
      out = dgamma(shift + (tinf-origin)) / normCons;
    }
   

  } else { // Asymptomatic infector and symptomatic infectors with missing symptom onset, they are infectious 2 days after their detection

    if ( infectionDate + 2.0 < studyPeriod && infectionDate + 2.0 < tinf) {
      out = dgamma(shift -3.0 + (tinf - infectionDate - 2.0)) / normCons;
    }
  }

  return out;
}


double cumulativeInfectivity(
                             double origin,
                             double infectionDate,
                             int infectionStatus,
                             double tinf,
                             int studyPeriod
                             )
{

  // Cumulative risk of infection
  // The probability density function was estimated by Aschcroft et al., 2020 (DOI: 10.4414/smw.2020.20336)
  double out = 0.0;
  double normCons = 1 - pgamma(shift - 3);

  if (infectionStatus == 1) { // Symptomatic infector with a date of symptom onset

    if (origin-3 < studyPeriod && origin-3 < tinf) { // The last data point for the recipient should be in ]origin-3, +inf[
      out = (pgamma(shift+(tinf-origin)) - pgamma(shift - 3.0)) / normCons;
    }

  } else { // Asymptomatic infector and symptomatic infectors with missing symptom onset, they are infectious 2 days after their detection

    if ( infectionDate + 2.0 < studyPeriod && infectionDate + 2.0 < tinf ){
      out = pgamma(shift - 3.0 +(tinf - infectionDate - 2.0)) - pgamma(shift - 3.0);
      out /= normCons;
    }
  }

  return out;
}

//----------------------------------------------------------------------
// Household class - Methods
//----------------------------------------------------------------------
Household::Household() : m_size(0), m_notInfected(0), m_startFollowUp(-1) {
    // Set size of vector attributes to 0
    m_indid.resize(0);
    m_onsetTime.resize(0.0);
    m_infected.resize(0);
    m_confCase.resize(0);
    m_studyPeriod.resize(0);
    m_vaccinationStatus.resize(0);
    m_age.resize(0);
    m_index.resize(0);
    m_isolation.resize(0);
    m_infTime.resize(0);
    m_cumLambda.resize(0);
    m_instLambda.resize(0);  
}


// Reinitialize household object
void Household::newHousehold() {
  m_size = 0;
  m_notInfected = 0;	
  m_startFollowUp = -1;
  m_indid.clear();
  m_onsetTime.clear();
  m_infected.clear();
  m_confCase.clear();
  m_studyPeriod.clear();
  m_vaccinationStatus.clear();
  m_age.clear();
  m_index.clear();
  m_isolation.clear();
  m_infTime.clear();
}

// Add an individual to the household
void Household::addIndividual(
  int indid, 
  double onsetTime, 
  int isCase, 
  int vaccinationStatus, 
  int age, 
  int studyPeriod, 
  int identifiedIndex, 
  int isolation
  ) {

	m_size += 1; // Update size of the household

	// Update vector attributes
	m_indid.push_back(indid);
	m_onsetTime.push_back(onsetTime);
	m_infected.push_back(isCase);
	m_studyPeriod.push_back(studyPeriod);
	m_vaccinationStatus.push_back(vaccinationStatus);
	m_age.push_back(age);
  	m_index.push_back(identifiedIndex); 
  	m_isolation.push_back(isolation);
	
	// Start of the follow-up of the household
	// It corresponds to the first symptom onset or the first detection date in the household
	if (m_startFollowUp < 0 ) m_startFollowUp = onsetTime;
	else m_startFollowUp = std::min(m_startFollowUp, onsetTime);

	// Add confirmed cases
	if (isCase > 0) { // 1: symptomatic cases, 2: asymptomatic cases, 3: symptomatic cases with unknown symptom onset
	    m_confCase.push_back(m_size-1);
	} else {
	    m_notInfected += 1;
	}

	// Initialize infection time
	m_infTime.push_back(1000.0);
}


// Display household
void Household::displayHH() {
  int ind;
  
  cout << "Identified index cases: ";
  for (ind=0; ind<m_size; ind++) cout << m_index[ind] << " "; 
  cout << endl;

  cout << "Onset Date: ";
  for (ind=0; ind<m_size; ind++) cout << m_onsetTime[ind] << " "; 
  cout << endl;

  cout << "Infection Date: ";
  for (ind=0; ind<m_size; ind++) cout << m_infTime[ind] << " "; 
  cout << endl;

  cout << "Study period: ";
  for (ind=0; ind<m_size; ind++) cout << m_studyPeriod[ind] << " "; 
  cout << endl;

  cout << "Vaccination status: ";
  for (ind=0; ind<m_size; ind++) cout << m_vaccinationStatus[ind] << " "; 
  cout << endl;

  cout << "Isolation behavior: ";
  for (ind=0; ind<m_size; ind++) cout << m_isolation[ind] << " "; 
  cout << endl;
 
  cout << "Infection status: ";
  for (ind=0; ind<m_size; ind++) cout << m_infected[ind] << " "; 
  cout << "\n\n";
}


// Access the infection status of a set of individuals 
std::vector<int> Household::getSpInfected(std::vector<int> index) const{
    
    unsigned i, k;
    std::vector<int> output;
    output.resize(0);

    for(i=0; i<index.size();i++) {
        k = index[i];
        output.push_back(m_infected[k]);
    }
    
    return output;
}


// Initialize the infection status of a single individual
void Household::initialInfTime(int index, std::mt19937_64& gen) {

    double incubPeriod(0.0);

    if (m_infected[index] == 1) {              // Symptomatic case with known symptom onset 
      // The incubation period was estimated by McAloon et al., 2020 (DOI: 10.1136/bmjopen-2020-039652)
      // It was truncated to [3-30] days so that all symptomatic individuals had the same infectivity profile
      incubPeriod = rlnorm(gen, mIncub, sdIncub); 
      while (incubPeriod < 3.0 || incubPeriod > 30.0) {
        incubPeriod = rlnorm(gen, mIncub, sdIncub);
      }
      m_infTime[index] = m_onsetTime[index] - incubPeriod;

    } else if (m_infected[index] > 1 ) {      // Asymptomatic case or symptomatic case with unknown symptom onset
      // Infection occurred in the 10 days before detection by PCR 
      m_infTime[index] = runif(gen, m_onsetTime[index] - maxPCRDetectability, m_onsetTime[index]);

    } else {		// Covid-free household members or household members with unknown final outcome
    	m_infTime[index] = 1000.0;
    }
}


// Change infection time of a single individual
void Household::setInfTime(int index, double infTime) {
	m_infTime[index] = infTime;
}


// At iteration 0 of the MCMC, compute the matrix of the pairwise instantaneous risk 
// of infection (m_instLambda) and the pairwise cumulative risk of infection (m_cumLambda) 
// between individuals in the household
// The first order of index corresponds to infectors and the second to infectees
void Household::compute_lambdas() {
  
  int infector, infectee;
  
  // Last day of follow-up for not infected contacts, symptom onset for symptomatic secondary cases, detection date by PCR for asymptomatic secondary cases
  double t1;  

  // Set m_cumLambda and m_instLambda size according to the household size
  m_cumLambda.resize(m_size);
  m_instLambda.resize(m_size);
  for (infector=0; infector<m_size; infector++) {
    m_cumLambda[infector].resize(m_size);
    m_instLambda[infector].resize(m_size);
  }

  // Compute cumulative and instantaneous person-to-person transmission rates
  for (infector=0; infector<m_size; infector++) {
    for (infectee=0; infectee<m_size; infectee++) {

      if ( m_onsetTime[infectee] == 1000.0 ) {
        t1 = m_studyPeriod[infectee];
      }else{
        t1 = m_infTime[infectee];
      }

      // Cumulative transmission rate from infector to secondary cases
      if (m_infTime[infector] != 1000.0 && m_infTime[infector] < t1) {

        m_cumLambda[infector][infectee] += cumulativeInfectivity(
          m_onsetTime[infector], 
          m_infTime[infector], 
          m_infected[infector], 
          t1, 
          m_studyPeriod[infectee]
        );

      }

      // Instantaneous transmission rate from infector to secondary cases
      if (m_infTime[infector] != 1000.0 && m_infTime[infectee] != 1000.0 && m_infTime[infector] < t1) {

        m_instLambda[infector][infectee] += infectivityProfile(
          m_onsetTime[infector], 
          m_infTime[infector], 
          m_infected[infector],
          m_infTime[infectee], 
          m_studyPeriod[infectee]
        );
      }

    }
  }

}

// Update m_cumLambda and m_instLambda matrices during data augmentation
void Household::update_lambdas(int ind) {
  
  int infector, infectee;
  double t1;

  // Update rows where ind is the infector
  for (infectee=0; infectee<m_size; infectee++) {
    // Reinitialize rows 
    m_cumLambda[ind][infectee] = 0;
    m_instLambda[ind][infectee] = 0;

    // Cumulative transmission rate for all household contacts
    if ( m_onsetTime[infectee] == 1000.0 ) {
      t1 = m_studyPeriod[infectee];
    }else{
      t1 = m_infTime[infectee];
    }

    if (m_infTime[ind] != 1000.0 && m_infTime[ind] < t1) {

      m_cumLambda[ind][infectee] = cumulativeInfectivity(
      m_onsetTime[ind], 
      m_infTime[ind], 
      m_infected[ind],
      t1, 
      m_studyPeriod[infectee]
      );
    }

    // Instantaneous transmission rate for secondary cases
    t1 = m_infTime[infectee];

    if (m_onsetTime[infectee] != 1000.0 && m_infTime[ind] != 1000.0 && m_infTime[ind] < t1) {
      m_instLambda[ind][infectee] = infectivityProfile(
      m_onsetTime[ind], 
      m_infTime[ind], 
      m_infected[ind], 
      t1, 
      m_studyPeriod[infectee]
      );
    }
  }


  // Update columns where ind is an infectee
  t1 = m_infTime[ind];

  for (infector=0; infector<m_size;infector++) {    
    // Reinitialize column 
    m_cumLambda[infector][ind] = 0;
    m_instLambda[infector][ind] = 0;

    // Cumulative transmission rate for all household contacts
    if (m_infTime[infector] != 1000.0 && m_infTime[infector] < t1) {

      m_cumLambda[infector][ind] = cumulativeInfectivity(
      m_onsetTime[infector], 
      m_infTime[infector], 
      m_infected[infector], 
      t1, 
      m_studyPeriod[ind]
      );
    }

    // Instantaneous transmission rate for secondary cases
    if (m_onsetTime[ind] != 1000.0 && m_infTime[infector] != 1000.0 && m_infTime[infector] < t1) {
      m_instLambda[infector][ind] = infectivityProfile(
      m_onsetTime[infector], 
      m_infTime[infector], 
      m_infected[infector], 
      t1, 
      m_studyPeriod[ind]
      );
    }
  }
}



// Sample a new infection time for a single individual
double Household::newInfTime(int index, std::mt19937_64& gen) {

    double infDate(0.0), incubPeriod(0.0);

    if (m_infected[index] == 1) {         	// Symptomatic case with known symptom onset 

      incubPeriod = rlnorm(gen, mIncub, sdIncub);
      while (incubPeriod < 3.0 || incubPeriod > 30.0) {
        incubPeriod = rlnorm(gen, mIncub, sdIncub);
      }
      infDate += m_onsetTime[index] - incubPeriod;

    } else if (m_infected[index] > 1) {	// Asymptomatic case or symptomatic case with unknown symptom onset
      infDate += runif(gen, m_onsetTime[index] - maxPCRDetectability, m_onsetTime[index]);

    } else {		// Covid-free household members or household members with unknown final outcome
    	infDate = 1000.0;
    }

    return infDate;
}

// Probability density of the incubation period for symptomatic cases and 
// the time from infection to detection for asymptomatic cases
double Household::pIncub(int index, double infTime) {

    double out = 0.0;

    if (m_infected[index] == 1) {       // Symptomatic case with known symptom onset 
      double normCons = plnorm(30.0) - plnorm(3.0);
      double incubPeriod = m_onsetTime[index] - infTime;

      if (incubPeriod >= 3.0) out += dlnorm(incubPeriod, mIncub, sdIncub) / normCons;

    } else {                      // Asymptomatic case or symptomatic case with unknown symptom onset
      out += 1/maxPCRDetectability;
    }

    return out;
}


// Log likelihood of the household
double Household::compute_log_lik(
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  double mainHHSize
  ) {

  // Identify the first case in the household
  int firstCase(0);
  auto it = std::min_element(m_infTime.begin(), m_infTime.end());
  double t0 = *it;
  for (int i =0; i < m_size; i++) {
    if ( m_infTime[i] == t0 ) firstCase = i;
  }

  // Log Likelihood
  double LL = 0.0;

  for (int i = 0; i < m_size; i++) { 
    
    if ( i == firstCase ) {                                           		// Contribution if the 1st case (incubation/detection period)
      LL += log( pIncub(i, m_infTime[i]) );

    } else if ( m_onsetTime[i] != 1000.0 && i != firstCase) {                 // Contribution of the secondary cases
        LL += log( S(i, t0, parameter, selectedParam, mainHHSize) );
      	LL += log( pInf(i, t0, parameter, selectedParam, mainHHSize) );
      	LL += log( pIncub(i, m_infTime[i]) );

    } else {									// Contribution of the contacts who were not infected
      LL += log( S(i, t0, parameter, selectedParam, mainHHSize) ); 

    }
  }
  
  return LL;
}


// Probability of infection at t1
double Household::pInf(
  int curr, 
  double t0, 
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  double mainHHSize
  ){

    // Instaneous risk of infection from community
    double alpha = parameter[0];

    // Instantaneous risk of infection from other household members
    double beta_i(0.0), beta(0.0);
    for (int ind=0; ind<m_size; ++ind) {
      
      beta_i = m_instLambda[ind][curr];
      
      // Relative infectivity of asymptomatic cases versus symptomatic cases
      if ( m_infected[ind] == 2) beta_i *= parameter[9];

      // Relative infectivity of vaccinated infectors
      if ( selectedParam[8] == 1 && m_vaccinationStatus[ind] == 1 ) beta_i *= parameter[8];

      beta += beta_i;
    }

    // Relative susceptibility
    if ( m_isolation[curr] ) {        // Infectees who isolated from the index case
	    if ( selectedParam[4] && m_age[curr] == 1 && m_vaccinationStatus[curr] ) beta *= parameter[4];
	    if ( selectedParam[5] && m_age[curr] == 1 && m_vaccinationStatus[curr] == 0 ) beta *= parameter[5];
	    if ( selectedParam[7] && m_age[curr] == 0 ) beta *= parameter[7];

    } else {                          // Infectees who did not isolate from the index case
    	if ( selectedParam[3] && m_age[curr] == 1 && m_vaccinationStatus[curr] ) beta *= parameter[3];
    	if ( selectedParam[6] && m_age[curr] == 0 ) beta *= parameter[6];

    }

    // Instantaneous per capita transmission rate 
    beta *= parameter[1];

    // If transmission hazard depends on household size
    if (selectedParam[2] == 1 && mainHHSize == 2.0) beta /= pow( m_size, parameter[2] );
    if (selectedParam[2] == 1 && mainHHSize != 2.0) beta /= pow( m_size/mainHHSize, parameter[2]);


    return 1-exp( -(alpha + beta) );
}




// Survival until from t0 to t1
double Household::S(
  int curr, 
  double t0, 
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  double mainHHSize
  ) {

  // If the individual is not infected, the infection time is coded as the  
  // last day of follow-up 
  double t1;
  if ( m_onsetTime[curr] == 1000.0 ) {
  	t1 = m_studyPeriod[curr];
  }else{
  	t1 = m_infTime[curr];
  }
  
  // Cumulative risk of infection from community
  double alpha = parameter[0] * (t1 - t0);

  // Cumulative risk of infection from other household members
  double beta(0.0), beta_i(0.0);
  for (int ind=0; ind<m_size; ++ind) {
  	beta_i = m_cumLambda[ind][curr];
	
	// Relative infectivity of asymptomatic cases versus symptomatic cases
	if ( m_infected[ind] == 2) beta_i *= parameter[9];
	
	// Relative infectivity of vaccinated infectors
	if ( selectedParam[8] == 1 && m_vaccinationStatus[ind] == 1 ) beta_i *= parameter[8];
	
	beta += beta_i;
   }

  // Relative susceptibility
  if ( m_isolation[curr] ) {        // Infectees who isolated from the index case
    if ( selectedParam[4] && m_age[curr] == 1 && m_vaccinationStatus[curr] ) beta *= parameter[4];
    if ( selectedParam[5] && m_age[curr] == 1 && m_vaccinationStatus[curr] == 0 ) beta *= parameter[5];
    if ( selectedParam[7] && m_age[curr] == 0 ) beta *= parameter[7];

  } else {                          // Infectees who did not isolate from the index case
    if ( selectedParam[3] && m_age[curr] == 1 && m_vaccinationStatus[curr] ) beta *= parameter[3];
    if ( selectedParam[6] && m_age[curr] == 0 ) beta *= parameter[6];

  }

  // Instantaneous per capita transmission rate
  beta *= parameter[1];
  
  // If transmission hazard depends on household size
  if (selectedParam[2] == 1 && mainHHSize == 2.0) beta /= pow(m_size, parameter[2]);
  if (selectedParam[2] == 1 && mainHHSize != 2.0) beta /= pow(m_size / mainHHSize, parameter[2]);
  
  return exp( -(alpha + beta ) );
}

