#include <math.h>
#include <float.h>

// experiment parameters
#define RANDOM_SEED 44635842
#define DURATION 400				// ms
#define N 200
#define NGROUPS 4
#define NUDGE_AT_US .02				// ms^-1
#define TRAININGCYCLES 300


//derived
#define TIMEBINS DURATION/DT

double *GE;							// unitless

void initDerivedParams() {
	GE = malloc(TIMEBINS * N * sizeof(double));
	for(int t = 0; t < TIMEBINS; t++) {
		for(int i = 0; i < N; i++ ) {
			if(i/(N/NGROUPS) <= t*DT/(DURATION/NGROUPS) && i/(N/NGROUPS) + 1 >= t*DT/(DURATION/NGROUPS)) *(GE + t*N + i) = NUDGE_AT_US * DT;
			else *(GE + t*N + i) = 0;
		}
	}	
}
