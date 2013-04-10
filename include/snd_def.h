
#ifndef ___INCLUDED_SND_DEF___
#define ___INCLUDED_SND_DEF___

#include <math.h>
#include <memory.h>
#include <sndfile.h>

#define def_vsize 256
#define def_sr 44100.f
#define def_cr def_sr/def_vsize
#define def_chans 1
#define	def_len 1024
#define pi 3.141592653

// Osc
float oscil(float amp, float freq, float* table, float* index, int len, float sr);
float osc(float *output, float amp, float freq, float *table, float *index,	int length, int vecsize, long sr);
float osci(float *output, float amp, float freq, float *table, float *index, float phase, int length, int vecsize, float sr);
float oscc(float *output, float amp, float freq, float *table, float *index, float phase, int length, int vecsize, float sr);

// Table generations
float* line_table(int brkpts, float* pts, int length);
float* exp_table(int brkpts, float* pts, int length);
void normalize_table(float* table, int length);

// Control signals
float line(float pos1, float dur, float pos2, int *cnt, float cr);
float expon(float pos1, float dur, float pos2, int *cnt, float cr);
float interp(float pos1, float dur, float pos2, double alpha, int *cnt, float cr);
float adsr(float maxamp, float dur, float at, float dt, float sus, float rt, int *cnt, float cr);

// Filters
float lowpass(float* sig, float freq, float *del, int vecsize, float sr);
float highpass(float* sig, float freq, float *del, int vecsize, float sr); 
float resonator(float* sig, float freq, float bw, float *del, int vecsize, float sr);
float bandpass(float *sig, float freq, float bw, float *del, int vecsize, int sr);
float balance(float *sig, float *cmp, float* del, float freq, int vecsize, float sr);
float delay(float *sig, float dtime, float *del, int *p, int vecsize, float sr);
float comb(float *sig, float dtime, float gain, float *delay, int *p, int vecsize, float sr);
float allpass(float *sig, float dtime, float gain, float *delay, int *p, int vecsize, float sr);
float vdelay(float *sig, float vdtime, float maxdel, float *delay, int *p, int vecsize, float sr);
float flanger(float *sig, float vdtime, float fdb, float maxdel, float *delay, int *p, int vecsize, float sr);


#endif