
#include "biquad.h"

BiQuad* BiQuad_create() {
	BiQuad* filter = (BiQuad*)malloc(sizeof(BiQuad));
	return filter;
}

void BiQuad_FlushDelays(BiQuad* filter) {
	filter->xz1 = 0;
	filter->xz2 = 0;
	filter->yz1 = 0;
	filter->yz2 = 0;
}

// Do the filter: given input xn, calculate output yn and return it
float BiQuad_do(BiQuad* filter, float xn) {
	// just do the difference equation: y(n) = a0x(n) + a1x(n-1) + a2x(n-2) - b1y(n-1) - b2y(n-2)
	float yn = filter->a0*xn + filter->a1*xz1 + filter->a2*xz2 - filter->b1*yz1 - filter->b2*yz2;

	// underflow check
	if(yn > 0.0 && yn < FLT_MIN_PLUS) yn = 0;
	if(yn < 0.0 && yn > FLT_MIN_MINUS) yn = 0;

	// shuffle the delays
	// Y delays
	filter->yz2 = filter->yz1;
	filter->yz1 = yn;

	// X delays
	filter->xz2 = filter->xz1;
	filter->xz1 = xn;

	return yn;
}