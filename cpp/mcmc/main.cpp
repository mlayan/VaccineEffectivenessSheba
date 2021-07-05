#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "Household.h"
#include "McmcObject.h"

using namespace std;

//----------------------------------------------------------------------
// Load data
//----------------------------------------------------------------------
std::vector<Household> buildData(std::string dataFile)
{

    // Variable with all the data
    std::vector<Household> output;
    output.resize(0);

    // Read file
    /*
    The data file should be a space-separated table with household ids
    sorted in ascending order and numbered from 0 to n-1 households
    */
    std::ifstream infile(dataFile.c_str());
    if(infile)
    {
        cout << "=============DATA===============" << endl;
        cout << "Read data file: " << dataFile << endl;

        int indid(0), lasthhid(0), currhhid(0), hhsize(0), isCase(0), vaccinationStatus(0), studyPeriod(0), age(0), identifiedIndex(0), isolation(0);
        double onsetTime = 0.0;
        int numberOfCase(0), numberOfHousehold(1), numberOfSubject(0), numberOfDay(0);
        //0: indid; 1: hhid; 2: hhsize; 3: onsetTime; 4: case; 5: test; 6: studyPeriod; 7: age; 8: identified index case; 9: isolation behavior from index case;

        Household currHH;

        for ( std::string line; std::getline(infile, line); )
        {
        	// Create a stringstream of the current line
        	std::istringstream in(line);
		
		// Store information in variables 
        	in >> indid >> currhhid >> hhsize >> onsetTime >> isCase >> vaccinationStatus >> studyPeriod >> age >> identifiedIndex >> isolation;

            // Append previous household to the list of households
            if (currhhid != lasthhid)
                {
                    numberOfHousehold++;
                    output.push_back(currHH);
                    currHH.newHousehold();
                    lasthhid = currhhid;
                }

            // Update household object and information parameters
            currHH.addIndividual(indid,
                                 onsetTime,
                                 isCase,
                                 vaccinationStatus,
                                 age,
                                 studyPeriod, 
                                 identifiedIndex, 
                                 isolation
                                 );
            numberOfSubject++;
            if (isCase != 0) numberOfCase++;
            numberOfDay = max(numberOfDay, studyPeriod);
        }

        // General information
        cout << "Number of households: " << numberOfHousehold << endl;
        cout << "Number of individuals: " << numberOfSubject << endl;
        cout << "Number of cases: " << numberOfCase<< endl;
        cout << "Number of days: " << numberOfDay << "\n\n";

        return(output);

    }else
    {
        cout << "ERROR: Cannot open the file." << endl;
        return(output);
    }
}


//----------------------------------------------------------------------
// Run mcmc
//----------------------------------------------------------------------
void runMCMC(McmcObject mcmc,
             std::string outputFile,
             int pas,
             std::vector<int> idOfSelectedParameter,
             std::vector<std::string> paramNames
)
{
	ofstream output(outputFile.c_str());

	int iteration, iter, parID;
	int numberOfIteration = int(mcmc.iteration() / pas);
	int nIterTimeInfection = mcmc.getNIterTimeInf();

	// Column names of the output file
	std::string colNames="iteration logLik ";
	for (size_t p=0; p < paramNames.size(); p++) colNames += paramNames[p] + " " + paramNames[p] + "_p " + paramNames[p] + "_a ";
	colNames += "data_p data_a";
	output << colNames << endl;

	// MCMC chain
	cout << "=============MCMC===============" << endl;

    	// Initial state
	mcmc.initialize_inf_time(); // Initialize infection time of all infected individuals
	mcmc.initial_log_lik(); 
	cout << "Initial log likelihood: " << mcmc.globalLogLik() << endl;
	
	output << "0 " << mcmc.globalLogLik() << " "; // Log likelihood
	for (size_t i = 0; i < mcmc.nParameters(); i++) {
	    output << mcmc.parameter(i) << " 0 0 ";
	    output << "0 0" << endl;
	}
	
	// Chain
	for (iteration = 0; iteration < numberOfIteration; iteration++)
	{
	    mcmc.resetMoves();

		for (iter = 0; iter < pas; iter++)
		{
			for (size_t selectedParameter = 0; selectedParameter < idOfSelectedParameter.size(); selectedParameter++)
			{
				parID = idOfSelectedParameter[selectedParameter];
				mcmc.update_parameter(parID);
			}

			// Data augmentation
			for (int i=0; i < nIterTimeInfection; i++) {
			    mcmc.update_augmented_data(); // Update log likelihood for all households
			}
		}

		// Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
		output << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood

		for (size_t i = 0; i < mcmc.nParameters(); i++) {
		    output << mcmc.parameter(i) << " " << mcmc.proposedMove(i) << " " << mcmc.acceptedMove(i) << " ";
		}
		
		output << mcmc.proposedMoveData() << " " << mcmc.acceptedMoveData() << endl;
	}
	
	output.close();
}



//----------------------------------------------------------------------
// Main function
//----------------------------------------------------------------------
int main(int argc, char **argv)
{

    // Arguments passed to main
    std::string model = argv[1];
    std::string hhSize = argv[2];
    int deltaParameter = std::stoi( argv[3] ); //0: default value to 1; 1: estimated 
    std::string chainID = argv[4];    
    double sdrInfVac = std::stod( argv[5] ); 
    double sdrS = std::stod( argv[6] ); 
    double asymp = std::stod( argv[7] ); //0.6, 0.95, 1 

    std::string vaccinationDefinition = argv[8]; // 1dose; 2doses
    std::string database = argv[9]; // name of the data base 

    double mainHHSize(2.0);
    if (argc >= 11) mainHHSize = std::stod( argv[10] ) ; // Main Household size 

    //==========Model parameters==========
    // Initial values
    std::vector<std::string> parameterNames = {"alpha", "beta", "delta", "rSAV", "rSAVI", "rSAI", "rSC", "rSCI", "rInfVac", "rAsymptomatic"};
    //                                            0        1       2       3        4       5     6        7        8 				9
    int numberOfParameters = parameterNames.size();
    std::vector<double> parameter(numberOfParameters);
    if (chainID == "1") parameter = {2e-4,  0.2,    0.5,    0.2,    0.6,    0.3,    1.5,    0.8,    0.3,	1.5	};
    if (chainID == "2") parameter = {4e-3,  0.5,    1.0,    0.3,    0.4,    0.9,    0.8,    1.0,    0.5,	0.4	};
    if (chainID == "3") parameter = {2e-5,  0.07,   1.5,    0.8,    0.53,   0.7,    1.0,    1.3,    1.0,	0.7	};

    // Parameters to infer according to model
    std::vector<int> selectedParameter(numberOfParameters, 1);

    if (model == "baseline") {
        selectedParameter[3] = 0;
        selectedParameter[4] = 0;
        selectedParameter[5] = 0;
        selectedParameter[6] = 0;
        selectedParameter[7] = 0;
        selectedParameter[8] = 0;
    }
    else if (model == "rS") {
        selectedParameter[8] = 0;
    } 
    else if (model == "rInfVac") {
    	selectedParameter[3] = 0;	
        selectedParameter[4] = 0;
        selectedParameter[5] = 0;
        selectedParameter[6] = 0;
        selectedParameter[7] = 0;
    }

    if (hhSize == "0") { 		// Power coefficient on the household size
    	selectedParameter[2] = 0; 
    } else {
    	if ( deltaParameter == 0 ) {
    		selectedParameter[2] = 0;
    		parameter[2] = 1; 
    	}
    }

    // Relative infectivity of asymptomatic cases
    selectedParameter[9] = 0; 
    parameter[9] = asymp;
    
    // Vector with the indices of the parameters to estimate
    int parameterNumber;
    std::vector<int> idOfSelectedParameter(0);
    for(parameterNumber=0; parameterNumber<numberOfParameters; parameterNumber++)
    {
        if(selectedParameter[parameterNumber]==1) idOfSelectedParameter.push_back(parameterNumber);
    }
	
	
   // Display model parametrization
    cout << "ID of selected parameters: "; 
    for (auto i = idOfSelectedParameter.begin(); i != idOfSelectedParameter.end(); ++i)
    	std::cout << parameterNames[*i] << ' ';
    cout << endl;

    cout << "Main household size: " << mainHHSize << "\n\n";


    //==========MCMC parameters==========
    size_t seed(20210329);
    int pas = 10;
    int numberOfIteration = 100000;
    int numberOfIterationTimeInfection = 1;

    // Standard deviation of the proposal distributions for the random walk
    std::vector<double> rateForRandomWalk(numberOfParameters);
    rateForRandomWalk = { 0.7, 0.6, 1.0, 2.2, 1.2, 1.3, 1.3, 1.8, 1.6, 0.5 };

    //==========Output files==========
    //Paths
    std::string fileName, pathData, pathOutput;
    pathData="/pasteur/sonic/homes/maylayan/MMMICovid/Israel/Data/";
    pathOutput="/pasteur/sonic/homes/maylayan/MMMICovid/Israel/Results/";

    //Data file
    std::string dataFile;
    std::string outputFile; 
    // File structure :     	0: indid; 1: hhid; 2: hhsize; 3: onsetTime; 4: case; 5: vaccinated; 6: studyPeriod; 7: age; 8: identified index case; 9: isolation behavior from index case;
    std::stringstream sssdrInfVac, sssdrS, ssAsymp;
    sssdrInfVac << std::fixed << std::setprecision(1) << sdrInfVac;
    sssdrS << std::fixed << std::setprecision(1) << sdrS;
    ssAsymp << std::fixed << std::setprecision(1) << asymp;
    dataFile=pathData + "2021_05_14_model_data_cpp_" + database + "_" + vaccinationDefinition + ".txt"; 
    outputFile=pathOutput + vaccinationDefinition + "/test3/mcmc_" + database + "_" + model + "_" + hhSize + "_" + std::to_string(int(maxPCRDetectability)) + "_" +  sssdrInfVac.str() + "_" + sssdrS.str() + "_" + ssAsymp.str() + "_" + chainID + ".txt";

    // Display names of input file and output file
    cout << "Input file: " << dataFile << endl;
    cout << "Output file: " << outputFile << "\n\n";

    //==========Build data==========
    // Load data
    std::vector<Household> hhData = buildData(dataFile);

    // Initialize MCMC object
    McmcObject mcmc(
        seed, 
        numberOfIteration, 
        hhData, 
        parameter, 
        selectedParameter, 
        rateForRandomWalk, 
        numberOfIterationTimeInfection, 
        mainHHSize, 
        sdrInfVac, 
        sdrS
    );

    //==========MCMC==========
    runMCMC(mcmc,
            outputFile,
            pas,
            idOfSelectedParameter,
            parameterNames
            );

    return 0;
}
