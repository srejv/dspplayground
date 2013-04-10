/** cross-synthesis main program

	usage: cross infile1 infile2 outfile dur(s)

	infile1: input1 (magnitudes) filename
	infile2: intput2 (phases) filename
	outfile: output filename
	dur: duration of input in seconds

	all files are supposed to be mono
*/

int main(int argc, char **argv)
{
	SNDFILE *fin, *fin2, *fout;
	SF_INFO input_info, input_info2;
	float *win, *in, *in2, *out, *spec, *spec2;
	int outsize, specsize, fftsize=1024, hopsize=256;
	int dataratio = fftsize/hopsize;
	int frames, extraframes;
	int read=0, written = 0, i, j, dur;
	float *buff;

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

	dur = input_info.frames > input_info2.frames ? input_info2.frames : input_info.frames;
	buff = (float*)malloc(sizeof(float)*100);
	win = (float*)malloc(sizeof(float)*fftsize);
	in = (float*)malloc(sizeof(float)*dur);
	in2 = (float*)malloc(sizeof(float)*dur);
	// number of full dataframes
	frames = ((dur*dataratio)/fftsize);
	// extra frames [not completely filled]
	extraframes = 1 + (frames*fftsize - dur*dataratio)/fftsize;
	// size of spectral data
	specsize = (frames + extraframes)*fftsize;
	spec = (float*)malloc(sizeof(float)*specsize);
	spec2 = (float*)malloc(sizeof(float)*specsize);
	outsize = specsize/dataratio;
	out = (float*)malloc(sizeof(float)*outsize);

	if(!(fout = sf_open(argv[3], SFM_WRITE, &input_info))) {
		printf("could not open %s \n", argv[3]);
		exit(-1);
	}

	// generera fönster
	for(i=0; i < fftsize; i++) win[i] = 0.5f - (float)(0.5*cos(i*twopi/fftsize));

	// läs in filer
	for(j=0; j < dur; j+=read) {
		read = sf_read_float(fin, buff, 100);
		for(i=0;i<read;i++) in[i+j] = buff[i];
	}
	for(j=0; j < dur; j+=read) {
		read = sf_read_float(fin2, buff, 100);
		for(i=0;i<read;i++) in2[i+j] = buff[i];
	}

	outsize = stft(in, win, spec, dur, fftsize, hopsize);
	stft(in2, win, spec2, dur, fftsize, hopsize);

	for(i=0; i < outsize; i+= fftsize)
		crosspec(&spec[i], &spec2[i], &spec[i], fftsize);

	dur = istft(spec, win, out, outsize, fftsize, hopsize);

	for(j=0; j < dur; j+= written) {
		for(i=0; i < 100 && j < dur; i++) {
			if(i+j < dur) buff[i] = out[i+j];
			else buff[i] = 0.f;
		}

		written = sf_write_float(fout, buff, i);
	}

	free(out);
	free(in);
	free(in2);
	free(spec);
	free(spec2);
	free(win);
	free(buff);

	sf_close(fout);
	sf_close(fin);
	sf_close(fin2);
	return 0;
}

void usage() {
	puts("\n\n usage: cross input1 input2 output \n");
}
