// envelope


#include <stdio.h>
#include <stdlib.h>
#include "snd_def.h"

/** simple synthesis program with envelopes
		
		Generates a sawtooth sound with frequency and amplitude
		envelopes. This program also shows the use of 
		interpolating oscillators.

		envelope sndfile.wav amp freq(Hz) dur(secs)
*/
int main(int argc, char const *argv[])
{
	SNDFILE *psf;
	float *buffer;
	float dur, amp, freq, *wave, ndx;
	int smps, cnt1=0, cnt2=0;

	if(argc == 5) {
		amp = (float)atof(argv[2]);
		freq = (float)atof(argv[3]);
		dur = (float)atof(argv[4]);
		smps =(int) (dur * def_cr);

		// allocate buffer & table memory
		buffer = new float[def_vsize];
		wave = sine_table();

		// open the file
		if(!(psf = soundout_open(argv[1]))) {
			printf("Error opening output file\n");
			exit(-1);
		}

		for(int i=0; i < smps; i++) {
			oscc(buffer,
				amp*adsr(1.f, dur, 0.05f, 0.1f, 0.7f, 0.2f, &cnt),
				expon(freq, dur/2, freq*2, &cnt2),
				wave, &ndx);
			sound.out(psf, buffer);
		}

		// close
		soundout_close(psf);
		delete[] buffer;
		delete[] wave;

		return 0;
	}
	else {
		printf("usage: envelope sndfile.wav amp freq(Hz) dur(s)\n");
		return 1;
	}
}