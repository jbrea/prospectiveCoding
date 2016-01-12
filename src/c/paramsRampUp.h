// Copyright 2016 Johanni Brea <johannibrea@gmail.com>
// 
// This file is part of prospectiveCoding.
// 
// prospectiveCoding is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
// 
// prospectiveCoding is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along with
// prospectiveCoding. If not, see http://www.gnu.org/licenses/.
// 
#include <math.h>
#include <gsl/gsl_rng.h>

// experiment parameters
#ifndef DURATION
	#define DURATION 2000		// timesteps in ms
#endif			
#ifndef NPRE	
	#define NPRE 2000			// number of presynaptic neurons
#endif
#define STIM_ONSET 1800			// onset of teaching stimulus in ms
#ifndef INPUT_AT_STIM				
	#define INPUT_AT_STIM .015	// input current in ms^-1
#endif
#ifndef TRAININGCYCLES
	#define TRAININGCYCLES 100
#endif
#define TAU_POSTS 10			// ms
#ifndef GAMMA_POSTS
	#define GAMMA_POSTS exp(- DT/TAU_POSTS)
#endif

#ifndef SEED_MAIN
	#define SEED_MAIN 230100138
#endif
#define SEED_PARAM 2389120

#define M_OU 1
#ifndef TAU_OU
	#define TAU_OU 400
#endif
#define N_OU 1
#define GAMMA_OU exp(- DT/TAU_OU)
#define S_OU sqrt(2 * DT / TAU_OU)

#define TIMEBINS DURATION/DT

double *GE, *GI, *PRE_ACT;

void initDerivedParams() {
	//gsl_rng *r;
	//gsl_rng_env_setup();
	//r = gsl_rng_alloc(gsl_rng_default);
	//gsl_rng_set(r, SEED_PARAM);
	 
	GE = malloc(TIMEBINS * sizeof(double)); 
	GI = malloc(TIMEBINS * sizeof(double)); 
	PRE_ACT = malloc(TIMEBINS * NPRE * sizeof(double)); 
	for(int t = 0; t < TIMEBINS; t++) {
		#if defined ROTATIONS
			double omega = 6.283185307179586/DURATION;
			GE[t] = .006 * (1 - sin(omega * t * DT) * sin(2 * omega * t * DT) * cos(4 * omega * t * DT)) * DT;
		#elif defined RAMPUPRATE
			if(t < 2*DURATION/DT/3) GE[t] = 0;
			else GE[t] = INPUT_AT_STIM * DT;
		#else 
			if(t > STIM_ONSET/DT) GE[t] = INPUT_AT_STIM * DT;
			else GE[t] = 0;
		#endif
				
		GI[t] = 0;
		#ifdef CURRENTPRED
			GI[t] = 4 * GE[t];
		#endif
		 
		#ifdef FROZEN_POISSON
			for(int i = 0; i < NPRE; i++) {
				if(gsl_rng_uniform(r) < .02 * DT) {
					PRE_ACT[NPRE * t + i] = 1;
				} else {
					PRE_ACT[NPRE * t + i] = 0;
				}
			}
		#elif defined RAMPUPRATE
			for(int i = 0; i < NPRE; i++) {
				if(t%((int)(DURATION/DT)/3) == 0) PRE_ACT[NPRE * t + i] = gsl_rng_uniform(r) * PHI_MAX/2 * DT;
				else PRE_ACT[NPRE * t + i] = PRE_ACT[NPRE * (t-1) + i];
			}
		#elif defined RAMPUPRATE_OU
			for(int i = 0; i < NPRE; i++) {
				if(t==0) PRE_ACT[i] = PHI_MAX/2 * DT*gsl_ran_gaussian_ziggurat(r,1);
				else PRE_ACT[t * NPRE + i] = runOU(PRE_ACT[(t-1) * NPRE + i], PHI_MAX/4 * DT, GAMMA_OU, PHI_MAX/2 * DT*S_OU);
			}
		#else 
			for(int i = 0; i < NPRE; i++) {
				if(i == t * DT) {
					PRE_ACT[NPRE * t + i] = 1;
				} else {
					PRE_ACT[NPRE * t + i] = 0;
				}
			}
		#endif
	}	
}
