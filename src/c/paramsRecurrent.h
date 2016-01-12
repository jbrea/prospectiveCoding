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
