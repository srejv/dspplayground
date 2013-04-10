
#include "mtapdelay.h"
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <math.h>

#define SAMPLING_RATE 44100.0

MTAP *new_mtap(long *tapindex, int ntaps, float maxtime, float *taptime, float *tapamp) {
	int i;
	MTAP *obj;	

	obj = (MTAP*)malloc(sizeof(MTAP) + (maxtime*SAMPLING_RATE * sizeof(float)));
	
	obj->ntaps = ntaps;
	obj->index = 0;
	obj->dl_length = maxtime*SAMPLING_RATE;
	obj->dl_space = 1.0f;
	obj->gain = 0.8;
	obj->fdb = 0.0;
	obj->wet = 0.4;

	for(i = 0; i < obj->ntaps; i++) {
		tapindex[i] = obj->dl_length - (long)( obj->dl_space * taptime[i] * SAMPLING_RATE + 0.5 );
	}

	for(i=0; i < obj->ntaps; i++) 
		tapamp[i] /= (float)sqrt((double)obj->ntaps);


	return obj;
}

/*
void mtap_update(MTAP *this, float *taptime, float *tapamp) {
	/* find starting position for each tap 
	int i;

	for(i=0; i < this->ntaps; i++) 
		this->tapamp /= sqrt((float)this->ntaps);

	for(i = 0; i < this->ntaps; i++) {
		this->tapindex[i] = this->dl_length - (long)( this->dl_space * this->taptime[i] * SAMPLING_RATE + 0.5 );
	}
}*/

void mtap_free(MTAP *this) {
	if(this) {

		if(this->delayline) {
			free(this->delayline);
		}

		free(this);
	}
}

void mtap_reset(MTAP *this) {
	this->index = 0;
	memset(this->delayline, 0, sizeof(float) * this->dl_length);
}

void mtap_process(MTAP *this, float *input, float *output, long blocksize, int chans, float wet, long *tapindex, float *taptime, float *tapamp) {
	int i, j;

	float dry = 1.0-wet;
	for(i = 0; i < blocksize; i++) {
		// in a processing loop
		for(j=0; j < this->ntaps; j++) {
			output[i] += this->delayline[tapindex[j]++] * tapamp[j];
			if(tapindex[j] == this->dl_length) {
				tapindex[j] = 0;
			}

			output[i] /= this->ntaps;
		}

		output[i] = wet*output[i] + input[i]*dry;

		// get new sample into delayline
		this->delayline[this->index++] = input[i];
		if(this->index == this->dl_length) {
			this->index = 0;
		}
	}
} 

