#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

gsl_rng *r;
#include "helper.h"
#include "paramsConditioning.h"

int main() {
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, RANDOM_SEED);
	initDerivedParams();
	
	//char filename[1024];
	//sprintf(filename, "data/conditioning_%d_%4f_%d.dat", TRAININGCYCLES, ETA, TAU_ALPHA);
	FILE *F = fopen(FILENAME, "w+b");
	
	gsl_vector *psp = gsl_vector_alloc(NPRE);
	gsl_vector *pspS = gsl_vector_alloc(NPRE);
	gsl_vector *sue = gsl_vector_alloc(NPRE);
	gsl_vector *sui = gsl_vector_alloc(NPRE);
	gsl_vector *pspTilde = gsl_vector_alloc(NPRE);
	gsl_vector *wB  = gsl_vector_alloc(NPRE);
	gsl_vector *wW  = gsl_vector_alloc(NPRE);
	gsl_vector *wR  = gsl_vector_alloc(NPRE);
	gsl_vector *ou  = gsl_vector_alloc(NPRE);
	gsl_vector *oum  = gsl_vector_alloc(NPRE);
	gsl_vector *pres  = gsl_vector_alloc(NPRE);
	double *pspP = gsl_vector_ptr(psp,0);
	double *pspSP = gsl_vector_ptr(pspS,0);
	double *sueP = gsl_vector_ptr(sue,0);
	double *suiP = gsl_vector_ptr(sui,0);
	double *pspTildeP = gsl_vector_ptr(pspTilde,0);
	double *wBP = gsl_vector_ptr(wB,0);
	double *wWP = gsl_vector_ptr(wW,0);
	double *wRP = gsl_vector_ptr(wR,0);
	double *ouP = gsl_vector_ptr(ou,0);
	double *oumP = gsl_vector_ptr(oum,0);
	double *presP = gsl_vector_ptr(pres,0);

	for(int i=0; i<NPRE; i++) {
		*(pspP+i) = 0;
		*(sueP+i) = 0;
		*(suiP+i) = 0;
		*(wBP+i) = gsl_ran_gaussian(r, .04) + .07;
		*(wWP+i) = gsl_ran_gaussian(r, .04) + .07;
		*(wRP+i) = gsl_ran_gaussian(r, .04) + .07;
	}
	
	double uB = 0, uVB = 0, rUB = 0, rVB = 0;
	double uW = 0, uVW = 0, rUW = 0, rVW = 0;
	double uR = 0, uVR = 0, rUR = 0, rVR = 0;
	
	
	int nOfRecordedPre = 50;
	int stateLength = 4 * nOfRecordedPre + 12;
	double *state[stateLength];
	for(int i = 0; i < nOfRecordedPre; i++) {
		*(state + 0*nOfRecordedPre + i) = wBP + i;
		*(state + 1*nOfRecordedPre + i) = wWP + i;
		*(state + 2*nOfRecordedPre + i) = wRP + i;
		*(state + 3*nOfRecordedPre + i) = presP + i;
	} 
	*(state + 4*nOfRecordedPre) = &uB;
	*(state + 4*nOfRecordedPre+1) = &uVB;
	*(state + 4*nOfRecordedPre+2) = &rUB;
	*(state + 4*nOfRecordedPre+3) = &rVB;
	*(state + 4*nOfRecordedPre+4) = &uW;
	*(state + 4*nOfRecordedPre+5) = &uVW;
	*(state + 4*nOfRecordedPre+6) = &rUW;
	*(state + 4*nOfRecordedPre+7) = &rVW;
	*(state + 4*nOfRecordedPre+8) = &uR;
	*(state + 4*nOfRecordedPre+9) = &uVR;
	*(state + 4*nOfRecordedPre+10) = &rUR;
	*(state + 4*nOfRecordedPre+11) = &rVR;
	
	
	double *IB, *IW, *IR = I1, *ou_t, uI;
	double IRf = 1;
	
	for( int s = 0; s < TRAININGCYCLES; s++) {
		if( (s <= TRAININGCYCLES - 9 && gsl_rng_uniform(r) > .5) || s%2==0) {
			ou_t = OU2; IB = I2; IW = I1; IRf = .5;
		} else {
			ou_t = OU1; IB = I1; IW = I2; IRf = 1; 
		}
		int interTrialInterval = gsl_rng_uniform_int(r, 7000 / DT);
		printf("%f\n", DT * interTrialInterval / 1000.);
		for( int t = 0; t < (US_END + 2000) / DT + interTrialInterval; t++) {
			for( int i = 0; i < NPRE; i++) {
				mixOUs(ouP + i, ou_t[t * NPRE + i], MIX[t], oumP + i);
				updatePre(sueP+i, suiP+i, pspP + i, pspSP + i, pspTildeP + i, *(presP + i) = spiking(DT * phi(*(oumP + i)), gsl_ran_flat(r,0,1))); 
			}
			updateMembrane(&uB, &uVB, &uI, wB, psp, IB[t], 0);
			updateMembrane(&uW, &uVW, &uI, wW, psp, IW[t], 0);
			updateMembrane(&uR, &uVR, &uI, wR, psp, IRf*IR[t], 0);
			//rUB = spiking(phi(uB), gsl_ran_flat(r,0,1)); rVB = phi(uVB);
			//rUW = spiking(phi(uW), gsl_ran_flat(r,0,1)); rVW = phi(uVW);
			//rUR = spiking(phi(uR), gsl_ran_flat(r,0,1)); rVR = phi(uVR);
			rUB = phi(uB); rVB = phi(uVB);
			rUW = phi(uW); rVW = phi(uVW);
			rUR = phi(uR); rVR = phi(uVR);
			for(int i = 0; i < NPRE; i++) {
				updateWeight(wBP + i, rUB, *(pspTildeP+i), rVB, *(pspSP+i));
				updateWeight(wWP + i, rUW, *(pspTildeP+i), rVW, *(pspSP+i));
				updateWeight(wRP + i, rUR, *(pspTildeP+i), rVR, *(pspSP+i));
			}
			if(s > TRAININGCYCLES - 9 && t < (US_END + 2000) / DT ) { 
				for(int i=0; i<stateLength; i++) {
					fwrite(*(state+i), sizeof(double), 1, F);
				}
			}
		}
	}
	
	gsl_vector_free(psp); gsl_vector_free(pspS); gsl_vector_free(wB); gsl_vector_free(wW); gsl_vector_free(wR);
	free(ou); free(oum); free(OU1); free(OU2); free(MIX); free(I1); free(I2);
	
	fclose(F); 
	
	return 0;
}
