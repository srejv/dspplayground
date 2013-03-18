
#include <stdlib.h>
#include <memory.h>
#include "delay.h"


DELAY * new_delay(float delay_time, float sample_rate, float gain, float feedback) {

	int length = (int)(delay_time * sample_rate);
	DELAY *d = (DELAY*)malloc(sizeof(DELAY) + length*sizeof(float));
	d->sr = sample_rate;
	d->delay = delay_time;
	d->wp = 0;
	d->length = length;
	d->gain = gain;
	d->feedback = feedback;
	memset(d->samples,0,length*sizeof(float));
	//d->samples = (float*)malloc(d->length*sizeof(float));
	
	return d;
}

void delay_free(DELAY** d) {
	if(d && *d && (*d)->samples){
		free((*d)->samples);
		free(*d);
		*d = NULL;
	}
}

void delay_tick(DELAY * d) {
	if(d->wp < d->length - 1) { 
		d->wp = d->wp + 1;
	} else {
		d->wp = 0;
	}
}

// Processera en frame (gör inline?)
float delay_processframe(DELAY * self, float in) {
	float out = self->samples[self->wp];
    self->samples[self->wp] = (self->gain * in) + (self->samples[self->wp] * self->feedback); 
    if(++(self->wp) >= self->length) self->wp = 0;

    return out;
}

// processera ett block
float delay_processblock(DELAY * self, float *in, float *out, unsigned long nFrames, int nChannels) {
	int i, j;
	for(i = 0; i < nFrames; i++) {
		for(j = 0; j < nChannels; j++) {
			*out++ = self->samples[self->wp];
			self->samples[self->wp] = (self->gain * *in++) + (self->samples[self->wp] * self->feedback);
			if(++(self->wp) >= self->length) self->wp = 0;
		}
	}

	return 0.0;
}
/*
float delay_vibrato_processframe(DELAY * self, float in) {

	float offset = self->wp - self->delay*self->sr;
	if(offset < 0.0) offset += self->delay*self->sr;
	float out = self->samples[offset];


}*/