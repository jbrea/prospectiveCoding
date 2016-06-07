rampUp:
	gcc  -std=c99 -O3 -D ETA=5 -D TRAININGCYCLES=2000 -D FILENAME_POST=\"data/rampUpPost.dat\" -D FILENAME_PRE=\"data/rampUpPre.dat\" -c src/c/rampUp.c -o rampUp.o
	gcc  rampUp.o -lgsl -lgslcblas -lm -o bin/rampUp
	rm rampUp.o
	time ./bin/rampUp
	
rampUpC:
	gcc  -std=c99 -O3 -D FILENAME_POST=\"data/rampUpPostC.dat\" -D FILENAME_PRE=\"data/rampUpPreC.dat\" -D CURRENTPRED -D ETA=50 -c src/c/rampUp.c -o rampUpC.o
	gcc  rampUpC.o -lgsl -lgslcblas -lm -o bin/rampUpC
	rm rampUpC.o
	time ./bin/rampUpC
	
rampUpFP:
	gcc  -std=c99 -O3 -D FROZEN_POISSON -D NPRE=500 -D TRAININGCYCLES=1000 -D RECORD_PREACT -D FILENAME_POST=\"data/rampUpFPPost.dat\" -D FILENAME_PRE=\"data/rampUpFPPre.dat\" -c src/c/rampUp.c -o rampUpFP.o
	gcc  rampUpFP.o -lgsl -lgslcblas -lm -o bin/rampUpFP
	rm rampUpFP.o
	time ./bin/rampUpFP
	
rampUpFPS:
	gcc  -std=c99 -O3 -D RECORD_PREACT -D TAU_OU=400 -D NPRE=500 -D SPIKING -D STOCHASTIC_US -D TRAININGCYCLES=1000 -D RECORD_PREACT -D FILENAME_POST=\"data/rampUpFPSPost.dat\" -D FILENAME_PRE=\"data/rampUpFPSPre.dat\" -D RAMPUPRATE_OU -D ETA=.2 -D SEED_MAIN=23411931 -c src/c/rampUp.c -o rampUpFP.o 
	gcc  rampUpFP.o -lgsl -lgslcblas -lm -o bin/rampUpFPS
	rm rampUpFP.o
	time ./bin/rampUpFPS
	
	
rampUpRate:
	gcc  -std=c99 -O3 -D RECORD_PREACT -D NPRE=500 -D SPIKING -D TRAININGCYCLES=1000 -D FILENAME_POST=\"data/rampUpRatePost.dat\" -D FILENAME_PRE=\"data/rampUpRatePre.dat\" -D INPUT_AT_STIM=.008 -D RAMPUPRATE -c src/c/rampUp.c -o rampUpRate.o
	gcc  rampUpRate.o -lgsl -lgslcblas -lm -o bin/rampUpRate
	rm rampUpRate.o
	time ./bin/rampUpRate
	
rampUpRateOU:
	gcc  -std=c99 -O3 -D RECORD_PREACT -D TAU_OU=400 -D NPRE=500 -D SPIKING -D TRAININGCYCLES=100 -D FILENAME_POST=\"data/rampUpRateOUPost.dat\" -D FILENAME_PRE=\"data/rampUpRateOUPre.dat\" -D RAMPUPRATE_OU -c src/c/rampUp.c -o rampUpRate.o
	gcc  rampUpRate.o -lgsl -lgslcblas -lm -o bin/rampUpRateOU
	rm rampUpRate.o
	time ./bin/rampUpRateOU

rotations:
	@for i in $$(seq 20 10 200); do \
		gcc -std=c99 -O3 -D ROTATIONS -D RECORD_PREACT -D TAU_ALPHA=$$i -D FILENAME_POST=\"data/rotationsPost_$$i.dat\" -D FILENAME_PRE=\"data/rotationsPre_$$i.dat\" -D TRAININGCYCLES=100 -D FROZEN_POISSON -c src/c/rampUp.c -o rotation$$i.o; gcc rotation$$i.o -lgsl -lgslcblas -lm -o bin/rotation$$i; rm rotation$$i.o; \
		time ./bin/rotation$$i; \
	done

delayedPairedAssociate:
	gcc  -std=c99 -O3 -D FILENAME=\"data/delayedPairedAssociate.dat\" -D TAU_ALPHA=2000 -c src/c/delayedPairedAssociate.c
	gcc  delayedPairedAssociate.o -lgsl -lgslcblas -lm -o bin/delayedPairedAssociate
	rm delayedPairedAssociate.o
	time ./bin/delayedPairedAssociate

recurrent:
	gcc  -std=c99 -O3 -D TAU_ALPHA=40 -c src/c/recurrent.c;
	gcc recurrent.o -lgsl -lgslcblas -lm -o bin/recurrent
	rm recurrent.o
	time ./bin/recurrent

makedirs:
	mkdir -p data bin

all: makedirs rampUp rampUpC rampUpFP rampUpRateOU rampUpFPS rampUpRate  rotations delayedPairedAssociate recurrent

