/*
	BiQuad struct
*/

#ifndef BIQUAD_H_
#define BIQUAD_H_

typedef struct BiQuad_t
{
	// Coefficients
	float a0, a1, a2, b1, b2;
	float c0, d0; // Not used atm (c0 = dry, d0 = wet)

	// Delays
	float xz1, xz2;
	float yz1, yz2;
} BiQuad;


BiQuad* BiQuad_create();
void BiQuad_FlushDelays(BiQuad* filter);
float BiQuad_do(BiQuad* filter, float xn);

#endif 