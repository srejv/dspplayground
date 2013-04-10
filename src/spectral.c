// Kassa versioner, vi ska använda FFTW tillfälligt. 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double twopi = 2*acos(-1.);

// dft takes an input signal *in of size N
// outputs sectrum *out with N pairs of 
// complex values [real, imag]
void dft(float *in, float*out, int N) {
	int i, k, n;
	for(i = 0; k = 0; k < N; i += 2, k++) {
		out[i] = out[i+1] = 0.f;
		for(int n = 0; n < N; n++) {
			out[i] += in[n]*cos(k*n*twopi/N);
			out[i+1] -= in[n]*sin(k*n*twopi/N);
		}
		out[i] /= N;
		out[i+1] /= N;
	}
}

// idft takes an input spectrum *in as N complex pairs
// outputs real signal *out, consisting of N samples
void idft(float *in, float *out, int N) {
	int n;
	for(n = 0; n < N; n++) {
		out[n] = 0.f;
		for( k = 0, i = 0; k < N; k++, i  += 2) {
			out[n] += in[i]*cos(k*n*twopi/N) - in[i+1] * sin( k*n*twopi/N);
		}
	}
}

void convol(float* inpulse, float* input, float* output,
	int impulse_size, int input_size) {

	float *impspec, *inspec, *outspec; 	// spectral vectors
	float *insig, *outsig, *overlap; 	// time domain vectors
	int fftsize = 1, convsize;			// transform and convolution sizes
	int overlap_size;					// overlap size
	int count, i, j; 					// counter and loop varibales

	overlap_size = impulse_size - 1;
	convsize = impulse_size + overlap_size;

	while(fftsize < convsize) fftsize *= 2;

	impspec =	(float*)malloc(sizeof(float) * fftsize);
	inspec =	(float*)malloc(sizeof(float) * fftsize);
	outspec =	(float*)malloc(sizeof(float) * fftsize);

	insig 	= (float*)malloc(sizeof(float) * fftsize);
	outsig	= (float*)malloc(sizeof(float) * fftsize);
	overlap = (float*)malloc(sizeof(float) * overlap_size);

	// get the impulse into the FFT input vector
	// pad with zeros
	for( i = 0; i < fftsize; i++) {
		if ( i < impulse_size ) insig[i] = impulse[i];
		else insig[i] = 0.f;
	}

	// Take the DFT of impulse
	fft(insig, impspec, fftsize);

	// processing loop
	for( i = count = 0; i < input_size+convsize; i++, count++) {

		// if an input block is ready
		if(count == impulse_size && i < (input_size+impulse_size)) {

			// copy overlapping block
			for( j = 0; j < overlap_size; j++ ) 
				overlap[j] = outsig[j+impulse_size];

			// pad input signal with zeros
			for( j = impulse_size; j < fftsize; j++)
				insig[j] = 0.f;

			// take the DFT of input signal block
			fft(insig, inspec, fftsize);

			// complex multiplication
			// first pair is re[0Hz] and re[Nyquist];
			outspec[0] = inspec[0]*impspec[0];
			outspec[1] = inspec[1]*impspec[1];

			// (a+ib) * (c+id) = (ac - bd) + (ad + bc)i
			for(j = 2; j < fftsize; j += 2) {
				outspec[j] = inspec[j]*impspec[j] - inspec[j+1]*impspec[j+1];
				outspec[j+1] = inspec[j]*impspec[j+1] + inspec[j+1]*impspec[j];
			}

			// IDFT of the spectral product
			ifft(otspec, outsig, fftsize);

			// zero the sample counter
			count = 0;
		}

		// get the input signal
		// stop when the input is finished
		if(i < input_size)
			insig[count] = input[i];

		// overlap-add output starts only 
		// after the first convolution operation
		if(i >= impulse_size)
			output[i-impulse_size] = outsig[ count ] +
				(count < overlap_size ? overlap[count] : 0);
	}

	free(overlap);
	free(outsig);
	free(insig);
	free(outspec);
	free(inspec);
	free(impspec);
}


#include "rfftw.h"
// these check if the FFTW okab was created
bool ft = true, ift = true;
static rfftw_plan forward, inverse; // these are the FFTW okabs

void fft(float *in, float *out, int N) {
	int i, k;
	if(ft) {
		// create the plan only once, first time it is called
		forward = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
		ft = false;
	}
	// copy input to the output
	for(i=0; i < N; i++) out[i] = in[i];
	// transform array out[] into in[]
	rfftw_one(forward, out, in);
	// move nyquist to pos 1
	out[0] = in[0]/N;
	out[1] = in[N/2]/N;
	// re-arrange array into [re, im] pairs
	for( i=2, k =1; i < N; k++, i+= 2) {
		out[i] = in[k]/N;
		out[i+1] = in[N-k]/N;
	}
}

void ifft(float *in, float *out, int N) {
	int i,k;
	if(ift) {
		// create the plan only once, first time it is called
		inverse = rfftw_create_plan(N, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
		ift = false;
	}
	// re-arrange array into the FFTW format
	for( i = 2, k = 1; i < N; k++, i+=2) {
		out[k] = in[i];
		out[N-k] = in[i+1];
	}
	// move Nyquist to pos N/2
	out[0] = in[0];
	out[N/2] = in[1];
	// transform array in[] into out[]
	rfftw_one(inverse, out, in);
	for(i=0; i < N; i++) out[i] = in[i];

}


int stft(float *input, float *window, float *output,
	int input_size, int fftsize, int hopsize) {

	int posin, posout, i;
	float *sigframe, *specframe;
	sigframe = (float*)malloc(sizeof(float)*fftsize);
	specframe = (float*)malloc(sizeof(float)*fftsize);

	for(posin = posout = 0; posin < input_size; posin += hopsize) {
		// window a signal frame
		for(i=0; i < fftsize; i++) {
			if(posin+i < input_size) 
				sigframe[i] = input[posin+i] * window[i];
			else sigframe[i] = 0;
		}
		// transform it
		fft(sigframe, specframe, fftsize);
		// output it
		for(i=0; i < fftsize; i++, posout++)
			output[posout] = specframe[i];
	}

	free(sigframe);
	free(specframe);
}

int stft(float* input, float* window, float* output, 
	int input_size, int fftsize, int hopsize) {

	int posin, posout, i, output_size;
	float *sigframe, *specframe;
	sigframe = (float*)malloc(sizeof(float)*fftsize);
	specframe = (float*)malloc(sizeof(float)*fftsize);
	output_size = input_size*hopsize/fftsize;

	for(posout=posin=0; posin < input_size; posin+=hopsize) {
		// window a signal frame
		for(i=0; i < fftsize; i++) {
			if(posin+i < input_size) {
				sigframe[i]  = input[posin+i]*window[i];
			}
			else {
				sigframe[i] = 0;
			}
		}

		// transform it
		fft(sigframe, specframe, fftsize);
		// output it
		for(i=0; i < fftsize; i++, posout++) 
			output[posout] = specframe[i];
	}

	free(sigframe);
	free(specframe);
	return posout;
}

int istft(float* input, float* window, float* output,
	int input_size, int fftsize, int hopsize) {

	int posin, posout, i, output_size;
	float *sigframe, *specframe;
	sigframe = (float*)malloc(sizeof(float)*fftsize);
	specframe = (float*)malloc(sizeof(float)*fftsize);
	output_size = input_size*hopsize/fftsize;

	for(posout=posin=0; posout < output_size; posout += hopsize) {
		// load in a spectral frame from input
		for(i=0; i < fftsize; i++, posin++) 
			specframe[i] = input[posin];

		// inverse-transform it
		ifft(specframe, sigframe, fftsize);

		//window it and overlap-add it
		for(i=0; i < fftsize; i++)
			if(posout+i < output_size)
				output[posout+i] += sigframe[i]*window[i];
	}

	free(sigframe);
	free(specframe);
	return output_size;
}

void crosspec(float *maginput, float *phasinput, float *output, int fftsize) {

	int i;
	float mag, phi;

	// take care of real-valued points at 0Hz and Nyquist
	output[0] = maginput[0];
	output[1] = maginput[1];
	for(i=2; i < fftsize; i+=2) {

		// get the magnitudes of one input
		mag = (float) sqrt(maginput[i]*maginput[i]+maginput[i+1]*maginput[i+1]);

		// get the phase of the other
		phi = (float)atan2(phasinput[i+1], phasinput[i]);

		// combine them and convert to rectangular form
		output[i] = (float) (mag*cos(phi));
		output[i+1] = (float) (mag*sin(phi));
	}
}

// Byt ut cos till sin så blir det en högpass! :D
void simplp(float *input, float *output, int fftsize) {

	int i,k;
	float mag, magin, phi;
	// The low-pass contour is 1 at 0Hz
	// and 0 at the Nyquist
	output[0] = 1.f;
	output[1] = 0.f;

	for( i=2 , k=1; i < fftsize; i+=2, k++) {

		// get the magnitueds of input
		magin = (float)sqrt(input[i]*input[i]+input[i+1]*input[i+1]);

		// apply the spectral contour
		mag = (float) cos(pi*k/fftsize)*magin;

		// get the phases
		phi = (float) atan2(input[i+1], input[i]);

		// convert to rectangular form
		output[i] = (float)(mag*cos(phi));
		output[i+1] = (float)(mag*sin(phi));
	}
}

void specreson(float *input, float *output, float scale,
	double angle, double radius, int fftsize) {

	int i, k;
	double sinw, cosw, cos2w, w, costheta, radsq, rad2;
	float re, im, div, rout, imout;

	costheta = cos(angle);
	radsq = radius * radius;
	rad2 = 2*radius;

	// 0 Hz and Nyquist taken care of
	output[0] = (float)((scale*input[0])/(1. - rad2*costheta + radsq));
	output[1] = (float)((scale*input[1])/(1. + rad2*costheta + radsq));

	for( i=2, k=1; i < fftsize; i+=2, k++ ) {

		w = (twopi*k)/fftsize;
		sinw = sin(w);
		cosw = cos(w);
		cos2w = cos(2*w);

		// real and imag parts of filter spectrum
		re = (float) (1. - rad2*costheta*cosw + radsq*cos2w);
		im = (float) (sinw*(rad2*costheta - 2*radsq*cosw));

		// complex division
		div = re*re + im*im;
		rout = (input[i]*re + input[i+1]*im)/div;
		imout = (input[i+1]*re - input[i]*im)/div;

		output[i] = scale*rout;
		output[i+1] = scale*imout;
	}
}

void specomb(float *input, float *output, float scale, float delay, double radius, int fftsize, float sr) {
	
	int i, k;
	double sinw, cosw, w, radsq, rad2;
	float re, im, div;
	radsq = radius*radius;
	rad2 = 2*radius;
	delay *= sr;

	// 0 Hz and Nyquist taken care of
	output[0] = (float)(input[0]*(1.-radius)/(1. - rad2 + radsq))*scale;
	output[1] = (float)(input[1]*-(1+radius)/(1. + rad2 + radsq))*scale;

	for( i=2, k=1; i < fftsize; i+=2, k++) {

		w = (delay*twopi*k)/fftsize;
		sinw = sin(w);
		cosw = cos(w);

		// real and imag parts of filter spectrum
		div = (float) (1. - rad2*cosw + radsq);
		re = (float) (cosw - radius) / div;
		im = (float) (sinw - rad2*cosw*sinw)/div;

		// complex multiplication
		output[i] = (input[i]*re - input[i+1]*im)*scale;
		output[i+1] = (input[i+1]*re + input[i]*im)*scale;
	}
}

