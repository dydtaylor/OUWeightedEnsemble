#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "twister.c"

#define RAND genrand_real1()
#define PI 3.14159265
#define ARRAYSIZE(a) (sizeof(a)/sizeof(a[0]))
#define NUM_ROWS(a) ARRAYSIZE(a)
#define NUM_COLS(a) ARRAYSIZE(a[0])

#define NBINSMAX 1000
#define BINCONTENTSMAXMAX 2000
#define ISIMMAXMAX 1000000 /* The largest this needs to be is NBINSMAX * BINCONTENTSMAXMAX,
	but we can probabilistically assume it is smaller */

struct paramsOUWeightedEnsemble{
	unsigned int tau; //In integer units of dt: e.g. if dt=.01, tau=.1, then this value should be .1/.01 = 10;
	unsigned int repsPerBin; // Target number of replicas per bin. Sometimes called "MTarg".
	unsigned int tauMax; //Number of weighted ensemble steps to do, max sim time.
			// If tauMax = 1000, tau = 10, dt = .005, then sim will end at t=tauMax*tau*dt = 50;
	unsigned int nBins;
	unsigned int fluxBin;
	double binDefs[NBINSMAX];
};

struct paramsOUDynamicsEngine{
	// Simulation parameters
	double dt; // Timestep of dynamics engine

	// Model parameters
	double tauSlow; // Timescale of OU process
	double sigmaX; // Standard deviation at steady-state of OU process
};

struct replicas{
	double sims[ISIMMAXMAX];
	double weights[ISIMMAXMAX];
	unsigned int binLocs[ISIMMAXMAX];
	unsigned int nBins;
	unsigned int binContentsMax[NBINSMAX];
	unsigned int binContents[BINCONTENTSMAXMAX][NBINSMAX];
	unsigned int iSimMax;
};

struct paramsOUWeightedEnsemble paramsWeOu;
struct paramsOUDynamicsEngine paramsDeOu;
struct replicas Reps;
FILE *errFile;

double dynamicsEngine(double x, struct paramsOUDynamicsEngine paramsDE, struct paramsOUWeightedEnsemble paramsWE){
	double U1, U2, Z1;
	double xOut;
	int nt;
	xOut = 0;
	for(nt = 0; nt < paramsWE.tau; nt++){ // JUN COMMENT: Should this go to tauMax? Or tau?
		U1 = RAND; // Box-Muller method for generating random normal variables
		U2 = RAND;
		Z1 = sqrt(-2*log(U1))*cos(2*PI*U2);
		xOut = x - 1/paramsDE.tauSlow * x * paramsDE.dt + paramsDE.sigmaX/sqrt(paramsDE.tauSlow) * sqrt(2*paramsDE.dt) * Z1;
	}
	
	if(isinf(xOut)){
				printf("NAN sim.\n");
	}
	return xOut;
}

int findBin(double rep, int nBins, struct paramsOUWeightedEnsemble paramsWE){
	/*Rewritten to only return a single integer for a single simulation (new binLocs for that index).
	In main, we will completely rewrite the table of contents each time.*/

	int binInside, iBin;
		for(iBin = 0; iBin < nBins; iBin++){
			if((rep <= paramsWE.binDefs[iBin+1]) && (rep > paramsWE.binDefs[iBin])){
				binInside = iBin;
				iBin = nBins;
			}
		}

	return binInside;
}

void initialDistOU(int nInit){
	/*
	This function currently sets the initial replicas to all be at x = 0.
	Obviously needs to be changed to actually draw random numbers if we want to take from the true dist.
	This most likely requires a loop?
	*/

	double startLocation = 0; //Can Call // JUN COMMENT ??
	int startBin, iRep, jBins;
	startBin = findBin(startLocation,paramsWeOu.nBins, paramsWeOu);

	for(jBins = 0; jBins<paramsWeOu.nBins; jBins++){
		Reps.binContentsMax[jBins] = 0;
	}

	for(iRep = 0; iRep < nInit; iRep++){
		Reps.sims[iRep] = startLocation;
		Reps.weights[iRep] = (double) 1/nInit;
		Reps.binLocs[iRep] = startBin;
		Reps.binContents[iRep][startBin] = iRep;
	}

	/*Want to rewrite this with the intention of modifying a struct without creating simsInit etc. Works once we malloc at the beginning of main.*/

	Reps.binContentsMax[startBin] = nInit;
	Reps.nBins = paramsWeOu.nBins;
	Reps.iSimMax = nInit - 1;
	return;
}

void splitMerge(){

	int binMin = 0;
	int binMax = Reps.nBins;
	unsigned int dummyInd;/*Just a dummy index that will hold the current index location in both split and merge loop*/
	int rowCol[3]; /* YET ANOTHER dummy that tells me which columns in the bincontents row will be combined together, with the third element giving the deleted column*/
	int mergeInd[2]; /* Array containing indices of 2 elements to be merged*/
	int keptInd[2]; /*Previous array, reordered s.t. the first index is the one that will be kept*/
	int splitInd; /* Index of replica to split*/
	int splitBin, mergeBin, repInBin, entryCheck1, entryCheck2, iSimMaxReplace; //Looping variables

	double p0; /*Doubles giving the merge probabilities*/
	double randPull; /*Double storing a pull from RAND*/
	//unsigned int binConPrev[BINCONTENTSMAXMAX][NBINSMAX];

	/*Merging loop*/
	for(mergeBin = binMin; mergeBin <binMax; mergeBin++){
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/

		while((Reps.binContentsMax[mergeBin]>paramsWeOu.repsPerBin) && Reps.binContentsMax[mergeBin] > 0){
			mergeInd[0] = Reps.binContents[0][mergeBin];
			mergeInd[1] = Reps.binContents[1][mergeBin];
			rowCol[0] = 0;
			rowCol[1] = 1;
			/*Find the locations of the two smallest weights and combine them together*/
			for(repInBin = 0; repInBin < Reps.binContentsMax[mergeBin];repInBin++){ //NUM_Rows needs to change to an appropriate BCM statement

				dummyInd = Reps.binContents[repInBin][mergeBin];
				/*If the weight of this index is greater than the weight in the first merge index,
				but smaller than the weight in the second index, then replace the first merge index with the new index

				Otherwise, if the weight of the dummy index is greater than the weight in the 2nd index, then replace the second
				index with the new index
				*/
				if((Reps.weights[dummyInd] < Reps.weights[mergeInd[0]])&&(Reps.weights[dummyInd] > Reps.weights[mergeInd[1]])){
					mergeInd[0] = dummyInd;
					rowCol[0] = repInBin;

				}
				else if((Reps.weights[dummyInd] < Reps.weights[mergeInd[1]])&& mergeInd[0] != dummyInd){
					mergeInd[1] = dummyInd;
					rowCol[1] = repInBin;
				}
			}

			if(mergeInd[0] == mergeInd[1]){
					fprintf(errFile, "Merge Error \n");
			}

			/*Decide which index to keep*/
			p0 = Reps.weights[mergeInd[0]] / (Reps.weights[mergeInd[0]]+Reps.weights[mergeInd[1]]);
			randPull = RAND;
			if(randPull<p0){ // JUN COMMENT: Do you need a variable randPull? Just use RAND, and, below, use else instead of else-if.
				keptInd[0] = mergeInd[0];
				keptInd[1] = mergeInd[1];
				rowCol[2] = rowCol[1];
			}
			else if(randPull > p0){
				keptInd[0] = mergeInd[1];
				keptInd[1] = mergeInd[0];
				rowCol[2] = rowCol[0];
			}

			/*Update weight of the kept index*/
			if( (isnan(Reps.weights[keptInd[0]])) || (isnan(Reps.weights[keptInd[1]]))){ // JUN COMMENT: Is it possible to use isnan()?
				fprintf(errFile, "WARNING: Moving NAN Weight \n");
			}
			Reps.weights[keptInd[0]] = Reps.weights[keptInd[0]] + Reps.weights[keptInd[1]];

			/*Replace the old simulation with the final non-NAN simulation*/

			if(isnan(Reps.weights[Reps.iSimMax])){
				fprintf(errFile, "WARNING: Moving NAN Weight \n");
			}

			int simMaxHolder;

			/*Find iSimMax in binContents and replace it with keptInd[1]*/
			/*For loop: iterates through the row of the bin ISM is in.
			If: Checks to see if the index of the row we're iterating through is ISM. If it is, it replaces it with KI[1], the deleted index.

			The reason we replace ISM with the deleted index is for the incrementing nature of keeping track of the sims.
			We aren't deleting ISM, we're moving it to the spot where the deleted sim was.

			This has problems when the deleted sim is ISM. When that happens, I will take the last element of the table of contents and
			move it to where ISM is in the table.*/

			for(iSimMaxReplace = 0; iSimMaxReplace < Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];iSimMaxReplace++){
				if(Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] == Reps.iSimMax){
					Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] = keptInd[1];
					simMaxHolder = iSimMaxReplace;
					iSimMaxReplace = Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];
				}
			}

			if (Reps.iSimMax ==keptInd[1]){
				Reps.binContents[simMaxHolder][Reps.binLocs[Reps.iSimMax]] = Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]] - 1;
			}

			Reps.sims[keptInd[1]] = Reps.sims[Reps.iSimMax];
			Reps.weights[keptInd[1]] = Reps.weights[Reps.iSimMax];
			Reps.binLocs[keptInd[1]] = Reps.binLocs[Reps.iSimMax];

			//Setting things to NAN is technically unnecessary
			/*Remove the duplicate non-NAN simulation at the end of the non-NANs*/
			Reps.sims[Reps.iSimMax] = NAN;
			Reps.weights[Reps.iSimMax] = NAN;
			Reps.binLocs[Reps.iSimMax] = NAN;
			Reps.iSimMax--;

			/*Reorganize the binContents matrix*/
			Reps.binContents[rowCol[2]][mergeBin] = Reps.binContents[Reps.binContentsMax[mergeBin] -1][mergeBin];
			Reps.binContents[-1+Reps.binContentsMax[mergeBin]][mergeBin] = NAN;
			Reps.binContentsMax[mergeBin]--;

			for(entryCheck1 = 0; entryCheck1 < Reps.binContentsMax[mergeBin]-1;entryCheck1++){
				for(entryCheck2 = entryCheck1 + 1; entryCheck2< Reps.binContentsMax[mergeBin]; entryCheck2++){
					if(Reps.binContents[entryCheck1][mergeBin]==Reps.binContents[entryCheck2][mergeBin]){
						fprintf(errFile, "ERROR: Duplicate entries in BC \n");
					}
				} // finished j loop through this bin's contents
			} // finished i loop through this bin's contents
		} // finished while-loop to check this bin
	} // finished merging loop through bins

	/*Splitting Loop*/
	for(splitBin = binMin; splitBin < binMax; splitBin++){
		while((Reps.binContentsMax[splitBin]<paramsWeOu.repsPerBin)&&(Reps.binContentsMax[splitBin]>0)){
			splitInd = Reps.binContents[0][splitBin];
			for(repInBin = 0; repInBin < Reps.binContentsMax[splitBin];repInBin++){
				dummyInd = Reps.binContents[repInBin][splitBin];
				if(Reps.weights[dummyInd]>Reps.weights[splitInd]){
					splitInd = dummyInd;
				}
			}
			Reps.sims[Reps.iSimMax+1] = Reps.sims[splitInd];
			Reps.weights[splitInd] = Reps.weights[splitInd] / 2;
			Reps.weights[Reps.iSimMax+1] = Reps.weights[splitInd];
			Reps.binLocs[Reps.iSimMax+1] = Reps.binLocs[splitInd];
			Reps.iSimMax++;
			if(Reps.iSimMax > ISIMMAXMAX){
				printf("ERROR: iSimMax out of bounds");
			}
			Reps.binContents[-1+Reps.binContentsMax[splitBin]][splitBin] = Reps.iSimMax;
			Reps.binContentsMax[splitBin]++;
		}
	} // finished splitting loop
	return;
}

double fluxes(){
	double fluxOut = 0;
	double newWeight, weightSum, weightPartialSum, scaleFactor;
	int jWeight, iReps, iSimMaxReplace, iSim;
	int nFlux = Reps.binContentsMax[paramsWeOu.fluxBin]; // Number of replicas in the fluxBin

	weightSum = 0;
	for(jWeight = 0; jWeight <= Reps.iSimMax; jWeight++){
		weightSum = weightSum + Reps.weights[jWeight];
	}

	// loop through replicas in flux bin and delete them
	for(iReps = 0; iReps < nFlux; iReps++){
		if(Reps.binLocs[iReps] == (paramsWeOu.fluxBin)){ // JUN COMMENT: Is this a validation check? Shouldn't all these be in the fluxbin?
			fluxOut += Reps.weights[Reps.binContents[iReps][paramsWeOu.fluxBin]]; // JUN COMMENT: Use += notation?
			newWeight = weightSum - fluxOut;
			scaleFactor = 1/newWeight;
			Reps.sims[Reps.binContents[iReps][paramsWeOu.fluxBin]] = Reps.sims[Reps.iSimMax];
			Reps.weights[Reps.binContents[iReps][paramsWeOu.fluxBin]] = Reps.weights[Reps.iSimMax];
			Reps.binLocs[Reps.binContents[iReps][paramsWeOu.fluxBin]] = Reps.binLocs[Reps.iSimMax];

			/*Find iSimMax in binContents and replace it with keptInd[1]*/
			for(iSimMaxReplace = 0; iSimMaxReplace < Reps.binContentsMax[Reps.binLocs[Reps.iSimMax]];iSimMaxReplace++){
				if(Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] == Reps.iSimMax){
					Reps.binContents[iSimMaxReplace][Reps.binLocs[Reps.iSimMax]] = Reps.binContents[iReps][paramsWeOu.fluxBin];
				}
			}

			Reps.sims[Reps.iSimMax] = NAN;
			Reps.weights[Reps.iSimMax] = NAN;
			Reps.binLocs[Reps.iSimMax] = NAN;
			Reps.iSimMax--;
		}
	}

	for(iReps = 0; iReps < Reps.binContentsMax[paramsWeOu.fluxBin]; iReps++){
		Reps.binContents[iReps][paramsWeOu.fluxBin] = NAN;
		weightSum = weightSum - Reps.weights[iReps];
	}
	Reps.binContentsMax[paramsWeOu.fluxBin] = 0;


	// JUN COMMENT: I don't understand why this while loop is necessary. Find total weight in one for-loop, then adjust everyone's weight in a second for-loop?
	if(fluxOut != 0){
		
			weightPartialSum = 0; //This variable was introduced because of the above while loop

			for(iSim = 0; iSim < Reps.iSimMax; iSim++){
				Reps.weights[iSim] = Reps.weights[iSim]*scaleFactor;
			}
			weightPartialSum = 0;
			for(jWeight = 0; jWeight <= Reps.iSimMax; jWeight++){
				weightPartialSum += weightPartialSum + Reps.weights[jWeight]; // JUN COMMENT: Use += notation?
			}
			weightSum = weightPartialSum;
			scaleFactor = 1/weightSum;
		//if(fabs(weightSum - 1)>1E-6){
			//printf("Weight Error: %E \n", fabs(weightSum - 1));
		//}
	}

	return fluxOut;
}

int main(int argc, char *argv[]){

	//Command line arguments: 1: File to record final simulations + weights
	//2: File to record fluxes
	//3: File to record errors
	//4: Bit for repeating old simulation. 1= repeat, 0 = new seed

	int tauQuarter, tauMax, rngBit, iBin, nWE;
	double binLoad;
	

	//Load parameters from files
	FILE *DEFile, *WEFile, *BINFile, *FLFile, *SIMFile;
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	BINFile = fopen("Bins.txt","r");
	errFile = fopen(argv[3], "w");

	fscanf(WEFile,"%i %i %i %i %i", &paramsWeOu.tau, &paramsWeOu.repsPerBin, &paramsWeOu.tauMax, &paramsWeOu.nBins, &paramsWeOu.fluxBin);
	paramsWeOu.fluxBin--;
	fscanf(DEFile, "%lf %lf %lf", &paramsDeOu.dt, &paramsDeOu.tauSlow, &paramsDeOu.sigmaX);

	for(iBin = 0; iBin<=(paramsWeOu.nBins);iBin++){ // JUN COMMENT: better counter name, like iBin.
		fscanf(BINFile, "%lf,",&binLoad);
		paramsWeOu.binDefs[iBin] = binLoad;
	}

	fclose(DEFile);
	fclose(WEFile);
	fclose(BINFile);

	printf("Parameters loaded\n");

	tauQuarter = paramsWeOu.tauMax / 4; //Sets number of steps for converging and data taking
	tauMax = paramsWeOu.tauMax;

	printf("Tau loops + flux vector made \n");

	rngBit = atoi(argv[4]);
	fprintf(errFile, "iseed=%lx\n", RanInitReturnIseed(rngBit));
    fclose(errFile);
	initialDistOU(paramsWeOu.repsPerBin);

	printf("Initial Distribution Made \n");
	printf("Initial Bin Location = %i \n", Reps.binLocs[0]);
	for(nWE = 0; nWE<tauQuarter; nWE++){
		printf("Tau Step: %i \n", nWE); //Show in stdout how far along in the program we are
		splitMerge();
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(int iSim = 0; iSim < Reps.iSimMax + 1;iSim++){
				if(isinf(Reps.sims[iSim])){
					printf("INF Error \n");
				}
				Reps.sims[iSim] = dynamicsEngine(Reps.sims[iSim],paramsDeOu, paramsWeOu);
				if(isinf(Reps.sims[iSim])){
					printf("INF Error \n");
				}
			
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim],Reps.nBins, paramsWeOu);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
	} // finished loop through taus for initial transient, e.g., first quarter

	for(nWE = tauQuarter; nWE<tauMax; nWE++){
		printf("Tau Step: %i \n", nWE);
		FLFile = fopen(argv[2],"a");
		fprintf(FLFile, "%E \n", fluxes());
		fclose(FLFile);
		splitMerge();
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(int iSim = 0; iSim < Reps.iSimMax;iSim++){
			Reps.sims[iSim] = dynamicsEngine(Reps.sims[iSim],paramsDeOu, paramsWeOu);
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim],Reps.nBins, paramsWeOu);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
	}
	 // finished loop through taus for non-transient (last 3 quarters)

	SIMFile = fopen(argv[1], "w");
	for(iBin = 0;iBin <= Reps.iSimMax; iBin++){ // JUN COMMENT: Better name for counter, like iBin.
		fprintf(SIMFile, "%E, %E \n", Reps.sims[iBin], Reps.weights[iBin]);
	}
	fclose(SIMFile);

	return 0;
}
