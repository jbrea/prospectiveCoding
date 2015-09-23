#include <math.h>
#include <float.h>

// experiment parameters
#define RANDOM_SEED 421251
#ifndef TRAININGCYCLES
	#define TRAININGCYCLES 200	// number of training cycles
#endif
#define DURATION 13000			// timesteps in ms
#define NPRE 2000				// number of presynaptic neurons
#define NUDGE_AT_US .005		// excitatory nudging conductance in ms^-1
#define TAU_OU 1000				// time constant of template OU in ms
#define TAU_OU2 100				// time constant of individual run OU in ms
#define CS_START 1000			//
#define CS_END 2000				// start and end of CS and US
#define US_START 3000			// in ms
#define US_END 4000				//
#define MIX_LOW1 0.1			// mixing parameter end of CS (unitless)
#define MIX_LOW2 0.2			// mixing parameter beginning of US (unitless)

//derived
double GAMMA_OU, GAMMA_OU2;
double SIGMA_OU, SIGMA_OU2;
int TIMEBINS;
double *OU1;
double *OU2;
double *MIX;
double *I1, *I2;

void mixOUs(double * restrict newOU, double oldOU, double mix, double * restrict mixedOU) {
	*mixedOU = runOU(*mixedOU, (1. - mix) * oldOU-.5 * mix, GAMMA_OU2, mix * SIGMA_OU2) ;
} 


void initDerivedParams() {
	 GAMMA_OU = exp(- DT/TAU_OU);
	 GAMMA_OU2 = exp(- DT/TAU_OU2);
	 SIGMA_OU = sqrt(2 * DT / TAU_OU);
	 SIGMA_OU2 = sqrt(2 * DT / TAU_OU2);
	 TIMEBINS = DURATION/DT;
	 
	 OU1 = malloc(TIMEBINS * NPRE * sizeof(double));
	 OU2 = malloc(TIMEBINS * NPRE * sizeof(double));	 
	 for( int i = 0; i < NPRE; i++) {
		OU1[i] = gsl_ran_gaussian_ziggurat(r,1);
		OU2[i] = gsl_ran_gaussian_ziggurat(r,1);
		for( int t = 1; t < TIMEBINS; t++) {
			OU1[t * NPRE + i] = runOU(OU1[(t-1) * NPRE + i], 0, GAMMA_OU, SIGMA_OU);
			OU2[t * NPRE + i] = runOU(OU2[(t-1) * NPRE + i], 0, GAMMA_OU, SIGMA_OU);
		}
	 }
	 
	 MIX = malloc(TIMEBINS * sizeof(double));
	 I1 = malloc(TIMEBINS * sizeof(double));
	 I2 = malloc(TIMEBINS * sizeof(double)); 
	 for(int t = 0; t < TIMEBINS; t++) {
		 I2[t] = 0;
		 if(t < US_START / DT || t > US_END / DT) I1[t] = 0; else I1[t] = NUDGE_AT_US * DT;
		 if(t < CS_START / DT || t > US_END / DT) {
			 MIX[t] = 1.;
		 } else if(t < CS_END / DT) {
			 MIX[t] = MIX[t-1] - DT * (1. - MIX_LOW1) / (CS_END - CS_START); 
		 } else if(t < US_START / DT) {
			 MIX[t] = MIX[t-1] + DT * (MIX_LOW2 - MIX_LOW1) / (US_START- CS_END);
		 } else MIX[t] = MIX[t-1] + DT * (1. - MIX_LOW2) / (US_END - US_START);
	 }	 
}
