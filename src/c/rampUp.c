#include <stdio.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

gsl_rng *r;
#include "helper.h"
#include "paramsRampUp.h"

#ifndef FILEPOST_FLAG
	#define FILEPOST_FLAG "wb"
#endif

int main() {
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, SEED_MAIN);
	initDerivedParams();
	
	FILE *postF = fopen(FILENAME_POST, FILEPOST_FLAG);
	FILE *preF = fopen(FILENAME_PRE, "wb");
	
	gsl_vector *psp = gsl_vector_alloc(NPRE);
	gsl_vector *pspS = gsl_vector_alloc(NPRE);
	gsl_vector *sue = gsl_vector_alloc(NPRE);
	gsl_vector *sui = gsl_vector_alloc(NPRE);
	gsl_vector *pspTilde = gsl_vector_alloc(NPRE);
	gsl_vector *w  = gsl_vector_alloc(NPRE);
	gsl_vector *pres  = gsl_vector_alloc(NPRE);
	#ifdef PREDICT_OU
		gsl_vector *ou = gsl_vector_alloc(N_OU);
		gsl_vector *preU = gsl_vector_calloc(NPRE);
		gsl_vector *wInput = gsl_vector_alloc(N_OU);
		gsl_matrix *wPre  = gsl_matrix_calloc(NPRE, N_OU);
		double *preUP = gsl_vector_ptr(preU,0);
		double *ouP = gsl_vector_ptr(ou,0);
		double *wInputP = gsl_vector_ptr(wInput,0);
		double *wPreP = gsl_matrix_ptr(wPre,0,0);
	#endif
	double *pspP = gsl_vector_ptr(psp,0);
	double *pspSP = gsl_vector_ptr(pspS,0);
	double *sueP = gsl_vector_ptr(sue,0);
	double *suiP = gsl_vector_ptr(sui,0);
	double *pspTildeP = gsl_vector_ptr(pspTilde,0);
	double *wP = gsl_vector_ptr(w,0);
	double *presP = gsl_vector_ptr(pres,0);

	for(int i=0; i<NPRE; i++) {
		*(pspP+i) = 0;
		*(sueP+i) = 0;
		*(suiP+i) = 0;
		#ifdef RANDI_WEIGHTS
			*(wP+i) = gsl_ran_gaussian(r, .1);
		#else
			*(wP+i) = 0;
		#endif
	}
	
	#ifdef PREDICT_OU
		for(int j=0; j < N_OU; j++) {
			*(ouP + j) = gsl_ran_gaussian(r, 1) + M_OU;
			*(wInputP + j) = gsl_ran_lognormal(r, 0., 2.)/N_OU/exp(2.)/2.;
			for(int i=0; i < NPRE; i++) *(wPreP + j*NPRE + i) = gsl_ran_lognormal(r, 0., 2.)/N_OU/exp(2.)/2.;
		}
	#endif
	
	double u = 0, uV = 0, rU = 0, rV = 0, uI = 0, rI = 0, uInput = 0;
	int on = 1;
	
	for( int s = 0; s < TRAININGCYCLES; s++) {
		#ifdef STOCHASTIC_US
			on = gsl_ran_bernoulli(r, .5);
		#endif
		for( int t = 0; t < TIMEBINS; t++) {
			#ifdef PREDICT_OU
				for(int i = 0; i < N_OU; i++) {
					*(ouP+i) = runOU(*(ouP+i), M_OU, GAMMA_OU, S_OU);
				}
				gsl_blas_dgemv(CblasNoTrans, 1., wPre, ou, 0., preU); 
			#endif
			for( int i = 0; i < NPRE; i++) {
				#ifdef SPIKING
					updatePre(sueP+i, suiP+i, pspP + i, pspSP + i, pspTildeP + i, *(presP + i) = spiking(PRE_ACT[t*NPRE + i], gsl_rng_uniform(r)));
				#elif defined PREDICT_OU
					//*(ouP+i) = runOU(*(ouP+i), M_OU, GAMMA_OU, S_OU);
					updatePre(sueP+i, suiP+i, pspP + i, pspSP + i, pspTildeP + i, *(presP + i) = DT * phi(*(preUP+i)));//spiking(DT * phi(*(preUP+i)), gsl_rng_uniform(r)));
				#else
					updatePre(sueP+i, suiP+i, pspP + i, pspSP + i, pspTildeP + i, *(presP + i) = PRE_ACT[t*NPRE + i]);
				#endif
			}
			#ifdef PREDICT_OU
				gsl_blas_ddot(wInput, ou, &uInput);
				GE[t] = DT * phi(uInput);
			#endif
			updateMembrane(&u, &uV, &uI, w, psp, on*GE[t], on*GI[t]);
			#ifdef POSTSPIKING
				rU = GAMMA_POSTS*rU + (1-GAMMA_POSTS)*spiking(DT * phi(u),  gsl_rng_uniform(r))/DT;
			#else
				rU = phi(u); 
			#endif
			rV = phi(uV); rI = phi(uI);
			for(int i = 0; i < NPRE; i++) {
				updateWeight(wP + i, rU, *(pspTildeP+i), rV, *(pspSP+i));
			}
			#ifdef TAUEFF
				if(s == TRAININGCYCLES - 1 && t < STIM_ONSET/DT) {
					fwrite(&rU, sizeof(double), 1, postF);
				}
			#else
				if(s%(TRAININGCYCLES/10)==0) {
					fwrite(&rU, sizeof(double), 1, postF);
					fwrite(GE+t, sizeof(double), 1, postF);
					fwrite(&rV, sizeof(double), 1, postF);
					fwrite(&rI, sizeof(double), 1, postF);
					fwrite(&u, sizeof(double), 1, postF);
				}
				if(s == TRAININGCYCLES - 1) {
					#ifdef RECORD_PREACT
						fwrite(PRE_ACT + t * NPRE, sizeof(double), 20, preF);
						//fwrite(ouP, sizeof(double), 20, preF);
						fwrite(presP, sizeof(double), 20, preF);
					#else
						fwrite(pspSP, sizeof(double), 1, preF);
						fwrite(pspTildeP, sizeof(double), 1, preF);
					#endif
				}
			#endif
		}
	}

	fclose(preF);
	fclose(postF);
	
	return 0;
}
