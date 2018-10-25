#include <stdio.h>
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

struct paramsOUWeightedEnsemble{
	int tau; 
	int repsPerBin;
	int tauMax;
	int nBins;
	int fluxBin;
	double binDefs[];
}

struct paramsOUDynamicsEngine{
	double dt;
	double tauSlow;
	double sigmaX;
}
/*
struct replicas{
	double sims[];
	double weights[];
	double binLocs[];
	int nBins;
	double binContentsMax[nBins];
	double binContents[][nBins];
	int iSimMax;
}*/

double dynamicsEngine(double x, struct paramsOUDynamicsEngine paramsDE, struct paramsOUWeightedEnsemble paramsWE){
	double U1, U2, Z1, Z2;
	double xOut;
	xOut = x;
	for(nt = 0; nt < tauMax; nt++;){
		U1 = RAND;
		U2 = RAND;
		Z1 = sqrt(-2*log(U1))*cos(2*PI*U2);
		Z2 = sqrt(-2*log(U1))*sin(2*PI*U2);
		xOut = xOut - 1/paramsDE.tauSlow * xOut * paramsDE.dt + paramsDE.sigmaX/sqrt(paramsDE.tauSlow) * sqrt(2*paramsDE.dt) * Z1;
	}
	return xOut;
}

void findBin(struct replicas currentReps, struct paramsOUWeightedEnsemble paramsWE){
	int alreadyIn;
	for(simIndex = 0; simIndex < currentReps.iSimMax; simIndex++){
		for(iBin = 0; iBin < nBins; iBin++){
			if((currentReps.sims[simIndex] < paramsWE.binDefs[iBin+1]) && (currentReps.sims[simIndex] > paramsWE.binDefs[iBin])){
				currentReps.binLocs[simIndex] = iBin;
				alreadyIn = 0;
				for(contentsLoop = 0; contentsLoop < currentReps.binContentsMax[iBin];contentsLoop++;){
					if(currentReps.binContents[contentsLoop][iBin] == simIndex){
						alreadyIn = 1;
					}
				}
				
				if(alreadyIn == 0){
					currentReps.binContentsMax[iBin]++;
					currentReps.binContents[currentReps.binContentsMax[iBin]][iBin] = simIndex;
				}
			}
		}
		
	}
	return;
}

void initialDistOU(struct paramsOUDynamicsEngine params, int nInit, int nMax, int nBins, double sims[], double weights[], double binLocs[], double binContentsMax[nBins],double binContents[][nBins], int iSimMax){
	/* 
	This function currently sets the initial replicas to all be at x = 0.
	Obviously needs to be changed to actually draw random numbers if we want to take from the true dist.
	This most likely requires a loop?
	*/
	double simsInit[nMax], weightsInit[nMax], binLocsInit[nMax], binContentsMaxInit[nBins];
	int binContsCols = 10*params.repsPerBin;
	double binContentsInit[binContsCols][nBins];
	double startLocation = 0;
	initialReplicas = malloc(sizeof(struct replicas) + 3*nMax*sizeof(double) + )
	int startBin;
	
	for(iBin = 0; iBin < nBins; iBin++){
		if(startLocation < params.binDefs[iBin+1] && startLocation > params.binDefs[iBin]){
			startBin = iBin;
		}
	}
	
	for(i = 0; i < nMax; i++){
		simsInit[i] = NAN;
		weightsInit[i] = NAN;
		binLocsInit[i] = NAN;
		for(k = 0;k < 10*params.repsPerBin;k++){
			binContentsInit[k][i] = NAN;
		}

	}
	
	for(jBins = 0; jBins<nBins; jBins++){
		for(iReps = 0; iReps < 10*params.repsPerBin; iReps++){
			binContentsInit[iReps][jBins] = NAN;
		}
		binContentsMaxInit[jBins] = 0;
	}
	
	for(j = 0; j < nInit; j++){
		simsInit[j] = startLocation;
		weightsInit[j] = 1/nInit;
		binLocsInit[j] = 0;
		binContentsInit[j][startBin] = j;
	}
	
	binContentsMaxInit[startBin] = nInit;
	sims = simsInit;
	weights = weightsInit;
	binLocs = binLocsInit;
	nBins = params.nBins;
	iSimMax = nInit;
	binContentsMax = binContentsMaxInit;
	binContents = binContentsInit;
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
	int loopIndex = 0; /* Dummy for checking how many times a merge / split loop operates */
	
	double p0, p1; /*Doubles giving the merge probabilities*/
	double randPull; /*Double storing a pull from RAND*/
	
	
	/* First, find the bins that we're looking through*/
	for(jBins = 0; jBins < nBins; jBins++){
		if((currentReps.BinContentsMax[jBins] > 0) &&(binMax > jBins) ){
			binMax = jBins;
		}
		if((currentReps.BinContentsMax[(nBins-jBins-1)] > 0) &&(binMin > (nBins-jBins-1))){
			binMin = nBins-jBins-1;
		}
	}
	if(binMax < binMin){
		fprintf(stdout, 'ERROR: binMax < binMin');
		return;
	}
	if(binMax == binMin){
		binMax++;
	}
	
	
	/*Merging loop*/
	for(mergeBin = binMin; mergeBin <binMax; mergeBin++){
		
		/*Initialize the merge indices with the first 2 indices stored in appropriate BinContents row*/
		mergeInd = {currentReps.BinContents[0][mergeBin],currentReps.BinContents[1][mergeBin]};
		
		
		while((currentReps.BinContentsMax[mergeBin]>params.repsPerBin) && currentReps.BinContentsMax[mergeBin] > 0){
			/*Find the locations of the two smallest weights and combine them together*/
			for(repInBin = 0; repInBin < NUM_ROWS(currentReps.BinContents);repInBin++){
				
				dummyInd = currentReps.BinContents[repInBin][mergeBin];
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
			for(iSimMaxReplace = 0; iSimMaxReplace < currentReps.BinContentsMax[currentReps.binLocs[currentReps.iSimMax]];iSimMaxReplace++;){
				if(currentReps.binContents[iSimMaxReplace][currentReps.binLocs[currentReps.iSimMax]] == currentReps.iSimMax){
					currentReps.binContents[iSimMaxReplace][currentReps.binLocs[currentReps.iSimMax]] = keptInd[1];
				}
			}
			
			/*Remove the duplicate non-NAN simulation at the end of the non-NANs*/
			currentReps.sims[currentReps.iSimMax] = NAN;
			currentReps.weights[currentReps.iSimMax] = NAN;
			currentReps.binLocs[currentReps.iSimMax] = NAN;
			currentReps.iSimMax--;
			
			/*Reorganize the binContents matrix*/
			currentReps.binContents[rowCol[2]][mergeBin] = currentReps.binContents[currentReps.binContentsMax[mergeBin]][mergeBin];
			currentReps.binContents[currentReps.binContentsMax[mergeBin]][mergeBin] = NAN;
			currentReps.binContentsMax--;
		}
	}
	
	
	/*Splitting Loop*/
	for(splitBin = binMin; splitBin < binMax; splitBin++){
		splitInd = currentReps.BinContents[0][splitBin];
		while((currentReps.BinContentsMax[splitBin]<params.repsPerBin)&&(currentReps.BinContentsMax[splitBin]>0)){
			for(repInBin = 0; repInBin < NUM_ROWS(currentReps.BinContents);repInBin++){
				dummyInd = currentReps.BinContents[repInBin][mergeBin];
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

double fluxes(struct replicas currentReps, struct paramsOUWeightedEnsemble paramsWE){
	double fluxOut = 0;
	
	
	
	for(iReps = 0; iReps < currentReps.iSimMax; iReps++){
		if(currentReps.binLocs[iReps] == paramsWE.fluxBin){
			fluxOut = fluxOut + currentReps.weights[iReps];
			currentReps.sims[iReps] = currentReps.sims[currentReps.iSimMax];
			currentReps.weights[iReps] = currentReps.weights[currentReps.iSimMax];
			currentReps.binLocs[iReps] = currentReps.binLocs[currentReps.iSimMax];
			
			/*Find iSimMax in binContents and replace it with keptInd[1]*/
			for(iSimMaxReplace = 0; iSimMaxReplace < currentReps.BinContentsMax[currentReps.binLocs[currentReps.iSimMax]];iSimMaxReplace++;){
				if(currentReps.binContents[iSimMaxReplace][currentReps.binLocs[currentReps.iSimMax]] == currentReps.iSimMax){
					currentReps.binContents[iSimMaxReplace][currentReps.binLocs[currentReps.iSimMax]] = iReps;
				}
			}
			
			currentReps.sims[currentReps.iSimMax] = NAN;
			currentReps.weights[currentReps.iSimMax] = NAN;
			currentReps.binLocs[currentReps.iSimMax] = NAN;
			currentReps.iSimMax--;
			
		}
	}
	
	
	for(iReps = 0; iReps < currentReps.binContentsMax[paramsWE.fluxBin]; iReps++;){
		currentReps.binContents[iReps][paramsWE.fluxBin] = NAN;
	}
	currentReps.binContentsMax[paramsWE.fluxBin] = 0;
	
	return fluxOut;
}

int main(){
	struct paramsOUWeightedEnsemble paramsWeOu;
	struct paramsOUDynamicsEngine paramsDeOu;
	
	struct replicas rep = initialDistOU(paramsDeOu, 500, 1000);
	double fluxes[paramsWeOu.tauMax];
	int tau1 = paramsWeOu.tauMax / 4;
	int tauM = paramsWeOu.tauMax;
	for(nWE = 1; nWE<tau1; nWE++){
		for(iSim = 0; iSim < rep.iSimMax;iSim++;){
			
		}
	}
}
