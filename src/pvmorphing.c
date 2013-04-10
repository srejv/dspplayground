#include <stdio.h>
#include "spec.h"
#include "fourier.h"

/** pv morphing main program



*/

int main(int argc, char **argv)
{
	SNDFILE *fin, *fin2, *fout;
	SF_INFO input_info, input_info2
	float *win, *in, *in2, *out, *spec, *spec2;	
	int outsize, specsize, fftsize=1024, hopsize = 256;
	int dataratio = fftsize/hopsize;
	int frames, extraframes;
	int read=0, written =0, i, j, dur;
	float *buff, sr;

	// Check arguments
	if(argc < 3) {
		printf("%s: unsufficient number of arguments(got %d, needed 2)", argv[0], argc-1);
		usage();
		exit(-1);
	}
	if(!(fin = sf_open(argv[1], SFM_READ, &input_info))) {
		printf("could not open %s", argv[1]);
		exit(-1);
	}
	if(!(fin2 = sf_open(argv[2], SFM_READ, &input_info2))) {
		printf("could not open %s", argv[2]);
		exit(-1);
	}
	if(!formats_equal(input_info, input_info2)) {
		sf_close(fin2);
		sf_close(fin);
		printf("%s and %s formats are not equal or files not mono\n", argv[1], argv[2]);
		exit(-1);
	}

	// load parameters
	dur = input_info.frames > input_info2.frames ? input_info2.frames : input_info.frames;
	sr = input_info.samplerate;
	buff = (float*)malloc(sizeof(float)*100);
	win = (float*)malloc(sizeof(float)*fftsize);
	in = (float*)malloc(sizeof(float)*dur);
	in2 = (float*)malloc(sizeof(float)*dur);
	// number of full dataframes
	frames = ((dur*dataratio)/fftsize);
	// extra frames (not completly filled)
	extraframes = 1 + (frames*fftsize - dur*dataratio)/fftsize;
	// size of spectral data
	specsize = (frames + extraframes) * fftsize;
	spec = (float*)malloc(sizeof(float)*specsize);
	spec2 = (float*)malloc(sizeof(float)*specsize);
	outsize = specsize / dataratio;
	out = (float*)malloc(sizeof(float)*outsize);

	if(!(fout = sf_open(argv[3], SFM_WRITE, &input_info))) {
		printf("could not open %s \n", argv[3]);
		exit(-1);
	}

	for(i=0; i < fftsize; i++) win[i] = 0.5f - (float)(0.5*cos(i*twopi/fftsize));

	for(j=0; j < dur; j+=read) {
		read = sf_read_float(fin, buff, 100);
		for(i=0; i < read; i++) in[i+j] = buff[i];
	}

	for(j=0; j < dur; j+= read) {
		read = sf_read_float(fin2, buff, 100);
		for(i=0; i < read; i++) in2[i+j] = buff[i];
	}

	outsize = pva(in, win, spec, dur, fftsize, hopsize, sr);
	pva(in2, win, spec2, dur, fftsize, hopsize, sr);

	for(i=0; i < specsize; i+= fftsize) {
		pvmorph(&spec[i], &spec2[i], &spec[i], (float)i/outsize, (float)i/outsize, 1024, (float)sr);
	}

	dur = psv(spec, win, out, specsize, fftsize, hopsize, sr);

	for(j=0; j < dur; j+= written) {
		for(i=0; i < 100 && j < dur; i++) {
			if(i+j < dur) buff[i] = out[i+j];
			else buff[i] = 0.f;
		}
		written = sf_write_float(fout, buff, i);
	}

	free(win);
	free(in);
	free(in2);
	free(out);
	free(spec);
	free(spec2);

	sf_close(fout);
	sf_close(fin);
	sf_close(fin2);
	return 0;
}

void usage() {
	puts("\n\n usage: morph input1 input2 output \n");
}
