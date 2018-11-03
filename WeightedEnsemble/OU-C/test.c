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
#define BINCONTENTSMAXMAX 1000
#define ISIMMAXMAX 1000000 /*NBINSMAX * BINCONTENTSMAXMAX */
#define TAUMAX 1000000

struct paramsOUWeightedEnsemble{
	int tau; //In units of dt: e.g. if dt=.01, tau=.1, then this value should be .1/.01 = 10;
	int repsPerBin;
	int tauMax; //Number of weighted ensemble steps to do, max sim time. If tauMax = 1000, tau = 10, dt = .005, then sim will end at t=tauMax*tau*dt = 50;
	int nBins;
	int fluxBin;
	double binDefs[];
};

struct paramsOUDynamicsEngine{
	double dt;
	double tauSlow;
	double sigmaX;
};

struct replicas{
	double sims[ISIMMAXMAX];
	double weights[ISIMMAXMAX];
	int binLocs[ISIMMAXMAX];
	int nBins;
	int binContentsMax[NBINSMAX];
	int binContents[BINCONTENTSMAXMAX][NBINSMAX];
	int iSimMax;
};

double dynamicsEngine(double x, struct paramsOUDynamicsEngine paramsDE, struct paramsOUWeightedEnsemble paramsWE){
	return x;
}

int findBin(double rep, int nBins, struct paramsOUWeightedEnsemble paramsWE){
	return nBins;
}

void initialDistOU(struct paramsOUWeightedEnsemble params, int nInit, struct replicas repsInit){
	/* 
	This function currently sets the initial replicas to all be at x = 0.
	Obviously needs to be changed to actually draw random numbers if we want to take from the true dist.
	This most likely requires a loop?
	*/
	
	double startLocation = 0; //Can Call 
	int startBin = findBin(startLocation,params.nBins,params); //

	
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


int main(int argc, char *argv[]){
	
	int userBins = atoi(argv[1]);
	if (userBins == 0){
		printf("Warning, no bins assigned");
	}
	struct paramsOUWeightedEnsemble *paramsWeOu = malloc(sizeof(*paramsWeOu) + userBins * sizeof(double));
	struct paramsOUDynamicsEngine paramsDeOu;
	struct replicas *Reps = malloc(sizeof(struct replicas));
	printf("%li", sizeof(*Reps));
	
	//struct replicas Reps;
	FILE *DEFile, *WEFile, *BINFile; //
	DEFile = fopen("dynamicsParams.txt","r");
	if(DEFile == NULL){
		printf("failDE");
	}
	WEFile = fopen("WEParams.txt","r");
	if(WEFile == NULL){
		printf("fail WE");
	}
	BINFile = fopen("Bins.txt","r");
	if(BINFile == NULL){
		printf("fail BIN");
	}
	fscanf(WEFile,"%i %i %i %i %i", &paramsWeOu->tau, &paramsWeOu->repsPerBin, &paramsWeOu->tauMax, &paramsWeOu->nBins, &paramsWeOu->fluxBin);
	fscanf(DEFile, "%lf %lf %lf", &paramsDeOu.dt, &paramsDeOu.tauSlow, &paramsDeOu.sigmaX);
	double binLoad[paramsWeOu->nBins+1];
	for(int j = 0; j<=paramsWeOu->nBins;j++){
		fscanf(BINFile, "%lf,",&binLoad[j]);
		paramsWeOu->binDefs[j] = binLoad[j];
		printf("%lf \n", paramsWeOu->binDefs[j]);
	}
	
	//int tau1 = paramsWeOu->tauMax / 4; //Sets number of steps for converging and data taking
	//int tauM = paramsWeOu->tauMax;
	double tauFluxes[TAUMAX] = {0};
	
	//initialDistOU(*paramsWeOu, paramsWeOu->repsPerBin,Reps); 
	/*
	for(int nWE = 0; nWE<tau1; nWE++){
		printf("Tau Step: %i \n", nWE); //Show in stdout how far along in the program we are
		splitMerge(Reps,*paramsWeOu);
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(int iSim = 0; iSim < Reps.iSimMax;iSim++){
			Reps.sims[iSim] = dynamicsEngine(Reps.sims[iSim],paramsDeOu, *paramsWeOu);
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim],Reps.nBins, *paramsWeOu);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
		
	}
	
	for(int nWE = tau1; nWE<tauM; nWE++){
		printf("Tau Step: %i \n", nWE); //Show in stdout how far along in the program we are
		splitMerge(Reps,*paramsWeOu);
		tauFluxes[tauM-nWE] = fluxes(Reps, *paramsWeOu);
		for(int iBin = 0; iBin < Reps.nBins; iBin++){
			Reps.binContentsMax[iBin] = 0;
		}
		for(int iSim = 0; iSim < Reps.iSimMax;iSim++){
			Reps.sims[iSim] = dynamicsEngine(Reps.sims[iSim],paramsDeOu, *paramsWeOu);
			Reps.binLocs[iSim] = findBin(Reps.sims[iSim],Reps.nBins, *paramsWeOu);
			Reps.binContents[Reps.binContentsMax[Reps.binLocs[iSim]]][Reps.binLocs[iSim]] = iSim;
			Reps.binContentsMax[Reps.binLocs[iSim]]++;
		}
	} */
	
	printf("Memory Sizes:  WEstruct %lu \n", sizeof(paramsWeOu));
	printf("Memory Sizes: DEstruct %lu \n", sizeof(paramsDeOu));
	printf("tau = %i \n",paramsWeOu->tau);
	printf("repsPerBin = %i \n",paramsWeOu->repsPerBin);
	printf("tauMax = %i \n",paramsWeOu->tauMax);
	printf("nBins = %i \n",paramsWeOu->nBins);
	printf("fluxBin = %i \n",paramsWeOu->fluxBin);
	
	printf("dt = %f \n", paramsDeOu.dt);
	printf("tauSlow = %f \n", paramsDeOu.tauSlow);
	printf("sigmaX = %f \n", paramsDeOu.sigmaX);
	
	fclose(DEFile);
	printf("close DE \n");
	fclose(WEFile);
	printf("clowe WE \n");
	fclose(BINFile);
	printf("close BIN \n");
	return 0;
}