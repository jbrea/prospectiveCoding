#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

gsl_rng *r;
#include "helper.h"
#include "paramsRecurrent.h"

int main() {
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, RANDOM_SEED);
	initDerivedParams();
	
	FILE *outF = fopen("data/recurrent.dat", "w+b");
	FILE *wDF = fopen("data/recurrentWD.dat", "w+b");
	FILE *wF = fopen("data/recurrentW.dat", "w+b");
	
	gsl_rng_set(r, 0);
	
	gsl_vector *psp = gsl_vector_calloc(N);
	gsl_vector *pspS = gsl_vector_calloc(N);
	gsl_vector *sue = gsl_vector_calloc(N);
	gsl_vector *sui = gsl_vector_calloc(N);
	gsl_vector *pspTilde = gsl_vector_calloc(N);
	gsl_matrix *w  = gsl_matrix_calloc(N, N);
	gsl_vector *u  = gsl_vector_calloc(N);
	gsl_vector *uV  = gsl_vector_calloc(N);
	gsl_vector *pre  = gsl_vector_calloc(N);
	gsl_vector *rU  = gsl_vector_calloc(N);
	gsl_vector *rV  = gsl_vector_calloc(N);
	double *pspP = gsl_vector_ptr(psp,0);
	double *pspSP = gsl_vector_ptr(pspS,0);
	double *sueP = gsl_vector_ptr(sue,0);
	double *suiP = gsl_vector_ptr(sui,0);
	double *pspTildeP = gsl_vector_ptr(pspTilde,0);
	double *wP = gsl_matrix_ptr(w,0,0);
	double *uP = gsl_vector_ptr(u,0);
	double *uVP = gsl_vector_ptr(uV,0);
	double *preP = gsl_vector_ptr(pre,0);
	double *rUP = gsl_vector_ptr(rU,0);
	double *rVP = gsl_vector_ptr(rV,0);

	
	int stateLength = 2 * N ;
	double *state[stateLength];
	for(int i = 0; i < N; i++) {
		*(state + 0*N + i) = preP + i;
		*(state + 1*N + i) = rUP + i;
	} 
	
	gsl_vector_view tmpv; double gE, wij, wd, uI;
	
	for( int s = 0; s < TRAININGCYCLES+1; s++) {
		wd = 0;
		for( int t = 0; t < TIMEBINS; t++) {
			for( int i = 0; i < N; i++) {
				updatePre(sueP+i, suiP+i, pspP + i, pspSP + i, pspTildeP + i, *(preP + i) = spiking(*(rUP+i), gsl_rng_uniform_pos(r))); 
				tmpv = gsl_matrix_row(w, i);
				if(s < TRAININGCYCLES - 1) gE = *(GE + t*N + i);
				else if(t > 7*TIMEBINS/8 && s < TRAININGCYCLES) gE = *(GE + i); //&& t < 2*TIMEBINS/3 + TIMEBINS/NGROUPS/5
				else gE = 0;
				updateMembrane(uP+i, uVP+i, &uI, &tmpv.vector, psp, gE, 0);
				*(rUP+i) = phi(*(uP+i)); *(rVP+i) = phi(*(uVP+i));
				for(int j = 0; j < N; j++) {
					wij = *(wP + i*N + j);
					if(i != j) {
						updateWeight(wP + i*N + j, *(rUP+i), *(pspTildeP+j), *(rVP+i), *(pspSP+j));
						wd += (wij -  *(wP + i*N + j)) * (wij -  *(wP + i*N + j));
					}
					if( s == TRAININGCYCLES - 1 && t == TIMEBINS - 1) fwrite(&wij, sizeof(double), 1, wF);
				}
			}
			if(s > TRAININGCYCLES - 4) {
				for(int i=0; i<stateLength; i++) fwrite(*(state+i), sizeof(double), 1, outF);
			}
		}
		fwrite(&wd, sizeof(double), 1, wDF);
	}

	fclose(outF); 
	fclose(wDF); 
	fclose(wF); 
	
	return 0;
}
