/*
	CDelay: implements a delay line of D samples with
	        output attenuator

	Can be used alone or as a base class.


	Copyright (c) 2010 Will Pirkle
	Free for academic use.

	Rewritten to C by Oscar Dragén
*/
#ifndef RDELAY_H_
#define RDELAY_H_

#include <math.h>

// underflow protection
#ifndef FLT_MIN_PLUS 
	#define FLT_MIN_PLUS          1.175494351e-38         /* min positive value */
#endif
#ifndef FLT_MIN_MINUS 
	#define FLT_MIN_MINUS        -1.175494351e-38         /* min negative value */
#endif

typedef struct RDelay_t {
	// member variables
	//
	// delay in samples; float supports fractional delay
	float fDelayInSamples;

	// output attenuation value, cooked
	float fOutputAttenuation;

	// read/write index values for circ buffer
	int nReadIndex;
	int nWriteIndex;

	// max length of buffer
	int nBufferSize;

	// sample rate (needed for other function)
	int nSampleRate;

	// delay in mSec, set by Parent Plug In
	float fDelay_ms;

	// Attenuation in dB, set by Parent Plug In
	float fOutputAttenuation_dB;

	// pointer to our circular buffer
	float* pBuffer;
} RDelay;

// constructor/destructor
RDelay * RDelay_create();
void RDelay_destroy(RDelay* this);  // base class MUST declare virtual destructor

// declare buffer and zero
void RDelay_init(RDelay* this, int nDelayLength);

// function to cook variables
void RDelay_cookVariables(RDelay* this);

// flush buffer, reset pointers to top
void RDelay_resetDelay(RDelay* this);

// set functions for Parent Plug In
void RDelay_setDelay_mSec(RDelay* this, float fmSec);
void RDelay_setSampleRate(RDelay* this, int nFs){ this->nSampleRate = nFs;};
void RDelay_setOutputAttenuation_dB(RDelay* this, float fAttendB);

// read the delay without writing or incrementing
float RDelay_readDelay(RDelay* this);

// read the delay at an arbitrary time without writing or incrementing
float RDelay_readDelayAt(RDelay* this, float fmSec);

// write the input and icrement pointers
void RDelay_writeDelayAndInc(RDelay* this, float fDelayInput);

// process audio
bool RDelay_processAudio(RDelay* this, float* pInput, float* pOutput);


#endif