#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>
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

        // while(std::getline(infile, line))
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

        // Add last household 
        output.push_back(currHH);

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
             /*int pasPrinting,*/
             std::vector<int> idOfSelectedParameter,
             std::vector<std::string> paramNames
)
{
	ofstream output(outputFile.c_str());

	int iteration, iter, parID;
	int numberOfIteration = int(mcmc.iteration() / pas);
	int nIterTimeInfection = mcmc.getNIterTimeInf();

	// Column names
	std::string colNames="iteration logLik ";
	for (size_t p=0; p < paramNames.size(); p++) colNames += paramNames[p] + " " + paramNames[p] + "_p " + paramNames[p] + "_a ";
	colNames += "data_p data_a";
	output << colNames << endl;

	// MCMC chain
	cout << "=============MCMC===============" << endl;

    // Initial state
    mcmc.initial_param_values(); // Initialize model parameters from prior
    cout << "Initial parameter values" << endl;
    for (size_t i = 0; i < idOfSelectedParameter.size(); i++) {
        cout << paramNames[idOfSelectedParameter[i]] << ": " << mcmc.parameter(idOfSelectedParameter[i]) << endl;
    }
	mcmc.initialize_inf_time(); // Initialize infection time of all infected individuals
	mcmc.initial_log_lik(); 
	cout << "Initial log likelihood: " << mcmc.globalLogLik() << endl;

    output << "0 " << mcmc.globalLogLik() << " "; // Log likelihood
    for (size_t i = 0; i < mcmc.nParameters(); i++)
        output << mcmc.parameter(i) << " 0 0 ";
    output << "0 0" << endl;

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

			// Augmented data
			for (int i=0; i < nIterTimeInfection; i++) {
			    mcmc.update_augmented_data(); // Update loglik at the household level
			}
		}

		// Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
		output << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood

		for (size_t i = 0; i < mcmc.nParameters(); i++)
            output << mcmc.parameter(i) << " " << mcmc.proposedMove(i) << " " << mcmc.acceptedMove(i) << " ";

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
    int hhSize = std::stoi(argv[2]);
    int deltaParameter = std::stoi( argv[3] ); //0: default value to 1; 1: estimated 
    std::string chainID = argv[4];

    double maxPCRDetectability = std::stod( argv[5] ); 
    if (maxPCRDetectability != 10.0 && maxPCRDetectability != 15.0 ) // Either 10 or 15
    	throw std::invalid_argument( "MaxPCRDetectability should be 10 or 15" );
    
    double sdrInfVac = std::stod( argv[6] ); 
    double sdrS = std::stod( argv[7] ); 
    double asymp = std::stod( argv[8] ); //0.6, 0.95, 1 

    std::string vaccinationDefinition = argv[9]; // 1dose; 2doses
    std::string runOnCluster = argv[10]; //0: desktop; 1: cluster
    std::string database = argv[11]; // name of the data base 

    double mainHHSize(2.0);
    if (argc >= 13) mainHHSize = std::stod( argv[12] ) ; // Main Household size in the case of French Guiana data  

    std::string fileName, infReduction;
    if (argc >= 14) {// For simulation studies, name of the input file and % of infectivity reduction of vaccinated cases 
    	fileName = argv[13]; 
    	infReduction = argv[14];
    }

    //==========Model parameters==========
    // Initial values
    std::vector<std::string> parameterNames = {"alpha", "beta", "delta", "rSAV", "rSAVI", "rSAI", "rSC", "rSCI", "rInfVac", "rAsymptomatic"};
    //                                            0        1       2       3        4       5     6        7        8 				9
    int numberOfParameters = parameterNames.size();
    std::vector<double> parameter(numberOfParameters);

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

    if (hhSize) { 			// Power coefficient on the household size
    	if ( deltaParameter == 0 ) {
            selectedParameter[2] = 0;
            parameter[2] = 1; 
        } 
    } else {
        selectedParameter[2] = 0;
    }

    // Relative infectivity of asymptomatic cases
    selectedParameter[9] = 0; 
    parameter[9] = asymp;
    
	int parameterNumber;
	std::vector<int> idOfSelectedParameter(0);
    for(parameterNumber=0; parameterNumber<numberOfParameters; parameterNumber++)
    {
        if(selectedParameter[parameterNumber]==1) idOfSelectedParameter.push_back(parameterNumber);
    }

    cout << "ID of selected parameters: "; 
    for (auto i = idOfSelectedParameter.begin(); i != idOfSelectedParameter.end(); ++i)
    	std::cout << parameterNames[*i] << ' ';
    cout << endl;

    cout << "Main household size: " << mainHHSize << "\n\n";


    //==========MCMC parameters==========
    size_t seed=0;
    switch (std::stoi(chainID)) {
        case 1:
            seed = 20210329;
            break;
        case 2:
            seed = 20210810;
            break;
        case 3:
            seed = 20210916;
            break;
        }
/*    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();*/
    cout << "Seed: " << seed << endl;

    int pas = 40;
    int numberOfIteration = 100000;
    int numberOfIterationTimeInfection = 1;

    // Variance of random walk
    std::vector<double> rateForRandomWalk(numberOfParameters);
    rateForRandomWalk = { 0.7,      0.6,    1.0,    2.2,    1.8,    1.3,    1.3,    1.8,    2.5,        0.5 };
    //                  "alpha", "beta", "delta", "rSAV", "rSAVI", "rSAI", "rSC", "rSCI", "rInfVac", "rAsymptomatic"    
    //==========Output files==========
    //Paths
    std::string pathData;
    std::string pathOutput;
    if (runOnCluster == "0")
    {
        pathData = "V:/maylayan/Israel/Data/";
        pathOutput = "V:/maylayan/Israel/Results/";
    } else
    {
        pathData="/pasteur/appa/homes/maylayan/MMMICovid/Israel/Data/";
        pathOutput="/pasteur/appa/homes/maylayan/MMMICovid/Israel/Results/";
    }

    //Data file
    std::string dataFile;
    std::string outputFile; 
    // File structure :     	0: indid; 1: hhid; 2: hhsize; 3: onsetTime; 4: case; 5: vaccinated; 6: studyPeriod; 7: age; 8: identified index case; 9: isolation behavior from index case;
    if (fileName.empty()) {
    	std::stringstream sssdrInfVac, sssdrS, ssAsymp;
    	sssdrInfVac << std::fixed << std::setprecision(1) << sdrInfVac;
    	sssdrS << std::fixed << std::setprecision(1) << sdrS;
        ssAsymp << std::fixed << std::setprecision(1) << asymp;
    	dataFile=pathData + "2021_05_14_" + database + "_" + vaccinationDefinition + ".txt";							//Input file 
		outputFile=pathOutput + vaccinationDefinition + "/mcmc_" + database + "_" + model + "_" + std::to_string(hhSize) + "_" + std::to_string(int(maxPCRDetectability)) + "_" +  
        sssdrInfVac.str() + "_" + sssdrS.str() + "_" + ssAsymp.str() + "_" + chainID + "_test.txt"; //Output file

    } else {
    	dataFile=pathData + vaccinationDefinition + "/simulation_" + infReduction + "/" + fileName;									//Input file 
		outputFile=pathOutput + vaccinationDefinition + "/simulation_" + infReduction + "/mcmc_" + chainID + "_" + fileName;			//Output file
    }
    
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
        sdrS,
        maxPCRDetectability
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
