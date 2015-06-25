#include <math.h>

// time step for forward Euler integration
#define DT .1								// in ms

// parameters of neural dynamics
#define GD 1.8								// in ms^-1 = mS/muF
#define GL .1								// in ms^-1 = mS/muF
#define GS .3								// in ms^-1 = mS/muF
#define EE 4.666666666666666666666667		// unitless
#define EI -.333333333333333333333333		// unitless
#ifndef TAU
	#define TAU 9							// in ms
#endif

// parameters of synaptic dynamics
#ifndef ETA
	#define ETA .5							// unitless
#endif
#ifndef TAU_ALPHA
	#define TAU_ALPHA 600					// ramp-up time constant in ms
#endif

// spiking parameters
#define PHI_MAX .06							// in kHz

//derived
#define	GD_DT GD * DT						// unitless
#define	GL_DT GL * DT						// unitless
#define	GAMMA_GL exp(- GL * DT)				// unitless
#define	GAMMA_GL_GD exp(- (GL+GD) * DT)		// unitless
#define	GAMMA_GS exp(- GS * DT)				// unitless
#ifdef CURRENTPRED
    #define GAMMA 0							// unitless
    #define ALPHA 1							// unitless
#else
	#define	GAMMA exp(- DT/TAU)
	#define	ALPHA (1. - exp(- DT * (1./TAU - 1./TAU_ALPHA)))
#endif
#define	C_PSP GL * GS/(GS-GL)				// unitless
#define ETA_DT ETA * DT						// ms

void updateMembrane(double * u, double * uV, double * uI, gsl_vector *w, gsl_vector *psp, double gE, double gI) {
	double v;
	gsl_blas_ddot(w, psp, &v); 
	*u = GAMMA_GL_GD * *u + GD_DT * v + gE * (EE - *u) + gI * (EI - *u);
	*uI = GAMMA_GL_GD * *uI +  gE * (EE - *u);
	*uV  = GAMMA_GL_GD * *uV + GD_DT * v;
}

void updatePre(double * sue, double * sui, double * psp, double * pspS, double * pspTilde, double pre) {
	*sue = GAMMA_GL * *sue + pre;
	*sui = GAMMA_GS * *sui + pre;
	*psp = C_PSP * (*sue - *sui);
	*pspS = *psp; //GAMMA_GL_GD * *pspS + GD_DT * *psp;
	*pspTilde = GAMMA * *pspTilde + *pspS;
}

void updateWeight(double * w, double rU, double pspTilde, double rV, double psp) {
	*w += ETA_DT * (ALPHA * rU * pspTilde - rV * psp);
}

double phi(double u) {
	if(u < 0) return 0;
	if(u > 1) return PHI_MAX;
	return PHI_MAX * u;
	//return PHI_MAX * DT / (1 + K * exp(BETA * (THETA - u)));
}

double spiking(double rate, double rn) {
	if( rn <=  rate) return 1.;
	return 0.;
}

double runOU(double ou, double m, double gamma, double s) {
	return gamma * ou + (1.-gamma) * m + s * gsl_ran_gaussian_ziggurat(r,1);
}
