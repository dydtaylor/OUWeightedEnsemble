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

#define NBINSMAX 100
#define BINCONTENTSMAXMAX 100
#define ISIMMAXMAX 10000 /*NBINSMAX * BINCONTENTSMAXMAX */
#define TAUMAX 100


struct paramsOUWeightedEnsemble{
	int tau; //In units of dt: e.g. if dt=.01, tau=.1, then this value should be .1/.01 = 10;
	int repsPerBin;
	int tauMax; //Number of weighted ensemble steps to do, max sim time. If tauMax = 1000, tau = 10, dt = .005, then sim will end at t=tauMax*tau*dt = 50;
	int nBins;
	int fluxBin;
	double binDefs[NBINSMAX];
};

struct paramsOUDynamicsEngine{
	double dt;
	double tauSlow;
	double sigmaX;
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

double dynamicsEngine(double x, struct paramsOUDynamicsEngine paramsDE, struct paramsOUWeightedEnsemble paramsWE){
	double U1, U2, Z1, Z2;
	double xOut;
	xOut = x;
	for(int nt = 0; nt < paramsWE.tauMax; nt++){
		U1 = RAND;
		U2 = RAND;
		Z1 = sqrt(-2*log(U1))*cos(2*PI*U2);
		Z2 = sqrt(-2*log(U1))*sin(2*PI*U2);
		xOut = xOut - 1/paramsDE.tauSlow * xOut * paramsDE.dt + paramsDE.sigmaX/sqrt(paramsDE.tauSlow) * sqrt(2*paramsDE.dt) * Z1;
	}
	return xOut;
}

int findBin(double rep, int nBins, struct paramsOUWeightedEnsemble paramsWE){
	/*Rewritten to only return a single integer for a single simulation (new binLocs for that index). 
	In main, we will completely rewrite the table of contents each time.*/
	int binInside;
		for(int iBin = 0; iBin < nBins; iBin++){
			printf("bin checked: %i \n", iBin+1);
			if((rep < paramsWE.binDefs[iBin+1]) && (rep > paramsWE.binDefs[iBin])){
				binInside = iBin+1;
				iBin = nBins;
				printf("Bin Found, bin = %i \n", binInside);
			}
		}
		
	return binInside;
}

void initialDistOU(struct paramsOUWeightedEnsemble params, int nInit, struct replicas repsInit){
	/* 
	This function currently sets the initial replicas to all be at x = 0.
	Obviously needs to be changed to actually draw random numbers if we want to take from the true dist.
	This most likely requires a loop?
	*/
	
	double startLocation = 0; //Can Call 
	int startBin;
	startBin = findBin(startLocation,params.nBins, params);

	
	/*117-119 is also unnecessary, sets values for locations that don't need to be set, 120 is necessary*/
	for(int jBins = 0; jBins<params.nBins; jBins++){
		repsInit.binContentsMax[jBins] = 0;
	}
	
	for(int j = 0; j < nInit; j++){
		repsInit.sims[j] = startLocation;
		repsInit.weights[j] = 1/nInit;
		repsInit.binLocs[j] = startBin;
		repsInit.binContents[j][startBin] = j;
	}
	
	/*Want to rewrite this with the intention of modifying a struct without creating simsInit etc. Works once we malloc at the beginning of main.*/
	
	repsInit.binContentsMax[startBin] = nInit;
	repsInit.nBins = params.nBins;
	repsInit.iSimMax = nInit;
	return;
}

void splitMerge(struct replicas currentReps, struct paramsOUWeightedEnsemble params){
	
	int binMax = 0;
	int binMin = currentReps.nBins;
	int dummyInd;/*Just a dummy index that will hold the current index location in both split and merge loop*/
	int rowCol[3]; /* YET ANOTHER dummy that tells me which columns in the bincontents row will be combined together, with the third element giving the deleted column*/
	int mergeInd[2]; /* Array containing indices of 2 elements to be merged*/
	int keptInd[2]; /*Previous array, reordered s.t. the first index is the one that will be kept*/
	int splitInd; /* Index of replica to split*/
	//int loopIndex = 0; /* Dummy for checking how many times a merge / split loop operates */
	
	double p0; /*Doubles giving the merge probabilities*/
	double randPull; /*Double storing a pull from RAND*/
	
	
	/*Merging loop*/
	for(int mergeBin = binMin; mergeBin <binMax; mergeBin++){
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/

		
		
		while((currentReps.binContentsMax[mergeBin]>params.repsPerBin) && currentReps.binContentsMax[mergeBin] > 0){
			mergeInd[0] = currentReps.binContents[0][mergeBin];
			mergeInd[1] = currentReps.binContents[1][mergeBin];
			/*Find the locations of the two smallest weights and combine them together*/
			for(int repInBin = 0; repInBin < currentReps.binContentsMax[mergeBin];repInBin++){ //NUM_Rows needs to change to an appropriate BCM statement
				
				dummyInd = currentReps.binContents[repInBin][mergeBin];
				/*If the weight of this index is greater than the weight in the first merge index,
				but smaller than the weight in the second index, then replace the first merge index with the new index
				
				Otherwise, if the weight of the dummy index is greater than the weight in the 2nd index, then replace the second 
				index with the new index
				*/
				if((currentReps.weights[dummyInd] < currentReps.weights[mergeInd[0]])&&(currentReps.weights[dummyInd] > currentReps.weights[mergeInd[1]])){
					mergeInd[0] = dummyInd;
					rowCol[0] = repInBin;
				
				}  
				else if(currentReps.weights[dummyInd] < currentReps.weights[mergeInd[1]]){
					mergeInd[1] = dummyInd;
					rowCol[1] = repInBin;
				
				}
			}
			
			/*Decide which index to keep*/
			p0 = currentReps.weights[mergeInd[0]] / (currentReps.weights[mergeInd[0]]+currentReps.weights[mergeInd[1]]);
			randPull = RAND;
			if(randPull<p0){
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
			currentReps.weights[keptInd[0]] = currentReps.weights[keptInd[0]] + currentReps.weights[keptInd[1]];
			
			/*Replace the old simulation with the final non-NAN simulation*/
			currentReps.sims[keptInd[1]] = currentReps.sims[currentReps.iSimMax];
			currentReps.weights[keptInd[1]] = currentReps.weights[currentReps.iSimMax];
			currentReps.binLocs[keptInd[1]] = currentReps.binLocs[currentReps.iSimMax];
			
			/*Find iSimMax in binContents and replace it with keptInd[1]*/
			for(int iSimMaxReplace = 0; iSimMaxReplace < currentReps.binContentsMax[currentReps.binLocs[currentReps.iSimMax]];iSimMaxReplace++){
				if(currentReps.binContents[iSimMaxReplace][currentReps.binLocs[currentReps.iSimMax]] == currentReps.iSimMax){
					currentReps.binContents[iSimMaxReplace][currentReps.binLocs[currentReps.iSimMax]] = keptInd[1];
				}
			}
			//Setting things to NAN is technically unnecessary
			/*Remove the duplicate non-NAN simulation at the end of the non-NANs*/
			currentReps.sims[currentReps.iSimMax] = NAN;
			currentReps.weights[currentReps.iSimMax] = NAN;
			currentReps.binLocs[currentReps.iSimMax] = NAN;
			currentReps.iSimMax--;
			
			/*Reorganize the binContents matrix*/
			currentReps.binContents[rowCol[2]][mergeBin] = currentReps.binContents[currentReps.binContentsMax[mergeBin]][mergeBin];
			currentReps.binContents[currentReps.binContentsMax[mergeBin]][mergeBin] = NAN;
			currentReps.binContentsMax[mergeBin]--;
		}
	}
	
	
	/*Splitting Loop*/
	for(int splitBin = binMin; splitBin < binMax; splitBin++){

		while((currentReps.binContentsMax[splitBin]<params.repsPerBin)&&(currentReps.binContentsMax[splitBin]>0)){
			splitInd = currentReps.binContents[0][splitBin];
			for(int repInBin = 0; repInBin < currentReps.binContentsMax[splitBin];repInBin++){
				dummyInd = currentReps.binContents[repInBin][splitBin];
				if(currentReps.weights[dummyInd]>currentReps.weights[splitInd]){
					splitInd = dummyInd;
				}
			}
			currentReps.sims[currentReps.iSimMax+1] = currentReps.sims[splitInd];
			currentReps.weights[splitInd] = currentReps.weights[splitInd] / 2;
			currentReps.weights[currentReps.iSimMax+1] = currentReps.weights[splitInd];
			currentReps.binLocs[currentReps.iSimMax+1] = currentReps.binLocs[splitInd];
			currentReps.iSimMax++;
			currentReps.binContents[currentReps.binContentsMax[splitBin]][splitBin] = currentReps.iSimMax;
			currentReps.binContentsMax[splitBin]++;
		}
	}
	return;
}

int main(int argc, char *argv[]){
	
	//int userBins = atoi(argv[1]); 
	
	struct paramsOUWeightedEnsemble paramsWeOu;
	
	struct paramsOUDynamicsEngine paramsDeOu;
	struct replicas Reps;
	int tau1, tauM, testInt;
	double binLoad;
	double tauFluxes[TAUMAX];
	
	//Load parameters from files
	FILE *DEFile, *WEFile, *BINFile; 
	DEFile = fopen("dynamicsParams.txt","r");
	WEFile = fopen("WEParams.txt","r");
	BINFile = fopen("Bins.txt","r");
	
	printf("Files opened \n");
	
	fscanf(WEFile,"%i %i %i %i %i", &paramsWeOu.tau, &paramsWeOu.repsPerBin, &paramsWeOu.tauMax, &paramsWeOu.nBins, &paramsWeOu.fluxBin);
	fscanf(DEFile, "%lf %lf %lf", &paramsDeOu.dt, &paramsDeOu.tauSlow, &paramsDeOu.sigmaX);

	for(int j = 0; j<=paramsWeOu.nBins;j++){
		fscanf(BINFile, "%lf,",&binLoad);
		paramsWeOu.binDefs[j] = binLoad;
	}
	fclose(DEFile);
	fclose(WEFile);
	fclose(BINFile);
	
	printf("Files loaded and closed \n");
	
	tau1 = paramsWeOu.tauMax / 4; //Sets number of steps for converging and data taking
	tauM = paramsWeOu.tauMax;
	
	printf("Tau loops + flux vector made \n");

	
	testInt = findBin(0.0, paramsWeOu.nBins, paramsWeOu);
	
	printf("Find Bin Location = %i \n",testInt);
	
	initialDistOU(paramsWeOu, paramsWeOu.repsPerBin, Reps);
	
	printf("Initial Distribution Made \n");
	return 0;
	//User specifies how many bins are going to be used, allows memory to be allocated appropriately given the FAM in WE struct
}
