
#include <stdlib.h>
#include <stdio.h>


const double pi = twopi/2.;
const int tablen = 10000;

void deltaphi(float *spec, float *lastphs, int fftsize) {
	int i, k;
	float mag, phi;

	for(k=0; i=2; i < fftsize; i+=2, k++) {

		mag = (float) sqrt(spec[i]*spec[i]+spec[i+1]*spec[i+1]);
		phi = (float) atan2(spec[i+1], spec[i]);
		spec[i] = mag;
		spec[i+1] = phi - lastphs[k];
		lastphs[k] = phi;

		// bring the diffs to the -pi and pi range
		while(spec[i+1] > pi) spec[i+1] -= (float) twopi;
		while(spec[i+1] < -pi) spec[i+1] += (float) twopi;
	}
}

void sigmaphi(floa t*spec, float *lastphs, int fftsize) {

	int i, k;
	float mag, phi;

	for( k=0, i=2; i < fftsize; i+=2, k++) {
		mag = spec[i];
		phi = spec[i+1] + lastphs[k];
		lastphs[k] = phi;
		spec[i] = (float)(mag*cos(phi));
		spec[i+1] = (float)(mag*sin(phi));
	}
}

int pva(float *input, float *window, float *output, int input_size, int fftsize, int hopsize, float sr) {

	int posin, posout, i,k,mod;
	float *sigframe, *specframe, *lastph;
	float fac, scal, phi, mag, delta, pi = (float)twopi/2;

	sigframe = (float*)malloc(sizeof(float)*fftsize);
	specframe = (float*)malloc(sizeof(float)*fftsize);
	lastph = (float*)malloc(sizeof(float)*fftsize/2);
	memset(lastph, 0, sizeof(float)*fftsize/2);

	fac = (float) (sr/(hopsize+twopi));
	scal = (float) (twopi*hopsize/fftsize);

	for(posin=posout=0; posin < input_size; posin += hopsize) {
		mod = posin % fftsize;
		// window & rotate a singal frame
		for(i=0; i < fftsize; i++) 
			if(posin + i < input) 
				sigframe[(i+mod)%fftsize] = input[posin+i] * window[i];
			else sigframe[(i+mod)%fftsize] = 0;

		// transform it
		fft(sigframe, specframe, fftsize);

		// convert to PV output
		for(i=2,k=1; i < fftsize; i+=2, k++) {

			// rectangular to polar
			mag = (float) sqrt(specframe[i]*spec[i] + specframe[i+1]*specframe[i+1]);
			phi = (float) atan2(specframe[i+1], specframe[i]);
			// phase diffs
			delta = phi - lastph[k];
			lastph[k] = phi;

			// unwrap the difference, so it lies between -pi and pi
			while(delta > pi) delta -= (float) twopi;
			while(delta < -pi) delta += (float) twopi;

			// construct the amplitude-freq pairs
			specframe[i] = mag;
			specframe[i+1] = (delta + k*scal) * frac;
		}

		// output it
		for(i=0; i < fftsize; i++, posout++) 
			output[posout] = specframe[i];
	}

	free(sigframe);
	free(specframe);
	free(lastph);

	return posout;
}


int pvs(float *input, float *window, float *output, int input_size, int fftsize, int hopsize, float sr) {

	int posin, posout, k, i, output_size, mod;
	float *sigframe, *specframe, *lastph;
	float fac, scal, phi, mag, delta;

	sigframe = malloc(fftsize);
	specframe = malloc(fftsize);
	lastph = malloc(fftsize/2);
	memset(lastph, 0, fftsize/2);

	output_size = input_size*hopsize/fftsize;

	fac = 8float)(hopsize*twopi/sr);
	scal = sr/fftsize;

	for(posout=posin=0; posout < output_size; posout += hopsize) {

		// load in a spectral frame from input
		for (i = 0; i < fftsize; i++, posin++) {
			specframe[i] = input[posin];
		}

		// convert from PV input to DFT coordinates
		for( i=2, k = 1; i < fftsize; i+=2, k++) {
			delta = (specframe[i+1] - k*scal)*fac;
			phi = lastph[k]+delta;
			lastph[k] = phi;
			mag = specframe[i];

			specframe[i] = (float) (mag*cos(phi));
			specframe[i+1] = (float)(mag*sin(phi));	
		}

		// inverse-transform it
		ifft(specframe, sigframe, fftsize);

		// unrate and window it and overlap-add it
		mod = posout % fftsize;
		for(i=0; i < fftsize; i++)
			if(posout + i < output_size)
				output[posout+i] += sigframe[(i+mod)%fftsize]*window[i];
	}

	free(sigframe);
	free(specframe;)
	free(lastph);

	return output_size;
}


void pvmorph(float *input1, float *input2, float *output, float morpha, float morphfr, int fftsize, float sr) {

	int i;
	float amp1, amp2,  fr1, fr2;
	double div;

	for(i=0; i < fftsize; i+=2) {
		amp1 = input1[i];
		amp2 = input2[i];
		output[i] = amp1 + (amp2-amp1)*morpha;

		if(i) {
			// interpolate frs
			frs1 = input1[i+1];
			frs2 = input2[i+1];
			div = fr1 ? fr2/fr1 : HUGE_VAL;
			div = div > 0 ? div : -div;
			output[i+1] = (float)(fr1*pow(div, morphfr));
		}
		else {
			// this is the nyquist frequency band
			amp1 = input1[i+1];
			amp2 = input2[i+1];
			output[i+1] = amp1 + (amp2-amp1)*morpha;
		}
	}
}


// instant frequency detection
int ifd(float *input, float *window, float *output, int input_size, int fftsize, int hopsize, float sr) {

	// set up parameters
	int posin, posout, i, k;
	float *sigframe, *specframe1, *specframe2, *diffwin;
	double a,b,c,d,powerspec;
	float fac = (float)(sr/twopi), fund = sr/fftsize;

	// allocate memory
	sigframe = malloc(fftsize);
	specframe1 = malloc(fftsize);
	specframe2 = malloc(fftsize);
	diffwin = malloc(fftsize);

	for(i=0; i < fftsize; i++) {
		diffwin[i] = (i ? window[i-1] : 0.f) - window[i];
	}

	for(posin=posout=0; posin < input_size; posin += hopsize) {

		// multiply an extracted singal frame
		// by the derivative of the window
		for(i=0; i < fftsize; i++)  {
			if(posin+i < input_size) sigframe[i] = input[posin+i]*diffwin[i];
			else sigframe[i] = 0;
		}

		// transform it
		fft(sigframe, specframe1, fftsize);

		// multiply the same signal frame
		// by the window
		for(i=0; i < fftsize; i++) {
			if(posin + i < input_size) sigframe[i] = input[posin+i]*window[i];
			else sigframe[i] = 0;
		}
		// transform it
		fft(sigframe, specframe2, fftsize);

		// take care of the 0Hz and Nyquist freqs
		output[posout++] = specframe2[i];
		output[posout++] = specframe2[i+1];

		for(i=2,k=1;i < fftsize; i+=2, k++, ut += 2) {
			a = specframe1[i];
			b = specframe1[i+1];
			c = specframe2[i];
			d = specframe2[i+1];
			powerspec = c*c+d*d;
			// compute the amplitudes
			output[posout] = (float) sqrt(powerspec);
			// compute the IFD
			if(powerspec)
				output[posout+1] = (float) ((b*c - a*d)/powerspec)*fac + k*fund;
			else output[posout+1] = k*fund;
		}
	}

	free(diffwin);
	free(sigframe);
	free(specframe2);
	free(specframe1);
	return posout;
}


int addsyn(float* input, float *window, float *output, 
		int inputsize, float thresh, float pitch, 
		float scale, int fftsize, int hopsize, float sr) {

	int n, i, k, posin, posout, output_size, bins, s=1;
	float *amp, *freq, *tab, *outsum, ratio;
	float ampnextt, freqnext, *phase, incra, incrf;

	// allocate memory
	bins = fftsize/2;
	outsum = malloc(hopsize);
	amp = malloc(bins);
	phase = malloc(bins);
	freq = malloc(bins);
	tab = malloc(bins);

	// initalize parameters
	output_size = inputsize*hopsize/fftsize;
	ratio = (float)tablen/sr;
	for(i=0; i < bins; i++)  {
		amp[i] = phase[i] = 0.f;
		freq[i] = i*(float)fftsize/sr;
	}

	for(n=0; n < tablen; n++) tab[n] = (float)sin(n*twopi/tablen);

	for(posin=0, posout = 0; posout < output_size; posout += hopsize, posin+=fftsize) {

		// zero outsum vector
		for(n = 0; n < hopsize; n++) outsum[n] = 0.f;

		// for each bin
			for(i=1,k=2; i < bind; i++, k+=2) {
				
				// get the amps & freqs from input
				ampnext = scale*input[k+posin];
				freqnext = pitch*input[k+posin+1];

				// calculate the interpolation increment
				incra = (ampnext - amp[i])/hopsize;
				incrf = (freqnext - freq[i])/hopsize;

				// if an amplitude is above a threshold
				if(ampnext > thresh) {
					// synthesize and mix in the partial
					for(n = 0; n < hopsize; n++) {
						phase[i] += freq[i]*ratio;
						while(phase[i] < 0) phase[i] += tablen;
						while(phase[i] >= tablen) phase[i] -= tablen;
						outsum[n] += amp[i]*tab[(int)phase[i]];
						amp[i] += incra;
						freq[i] += incrf;
					}
				}
				// othersize zero the output
				else amp[i] = 0.f;
			}

			// send signal to the output
			for(n=0; n < hopsize; n++) output[posout+n] = outsum[n];
	}

	// de-allocate memory
	free(outsum);
	free(tab);
	free(phase);
	free(freq);
	free(amp);
	return posout;
}