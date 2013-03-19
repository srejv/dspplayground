
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


// osc

float oscil(float amp, float freq, float* table, float* index, int len, float sr) {
	float out;
	out = amp*table[(int) *index];
	*index += freq*len/sr;
	while(index >= len) index -= len;
	while(index < 0) index += len;
	return out;
}

/** truncating table-lookup oscillator.

	output: output signal vector
	amp: amplitude
	freq: frequency
	table: function table
	index: lookup index
	length: function talbe length
	vecsize: signal vector size
	sr: sampling rate

	returns: first output sample
*/
float osc(float *output, float amp, float freq, 
	float *table, float *index,
	int length, int vecsize, long sr) {
	// increment
	float incr = freq*length/sr;

	// processing loop
	for(int i=0; i < vecsize; i++) {
		// truncated lookup
		output[i] = amp*table[(int)(*index)];
		*index += incr;
		while(*index >= length) *index -= length;
		while(*index < 0) *index += length;
	}

	return *output;
}

// osc usage
for(int i=0; i <dur; i++) {
	osc(buffer, amp, freq, wave, &ndx);
	soundout(psf, buffer);
}

// modulate the frequency. WOAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA OMG OMG OMG ;DDDD
for(int i=0; i < dur; i++) {
	osc(buffer, amp,
		freq+osc(&cs, 10.f, 5.f, wave, &ndx2, def_len. 1, def_cr),
		wave, &ndx);
	soundout(psf, buffer);
}

float osci(float *output, float amp, float freq, 
	float *table, float *index, float phase,
	int length, int vecsize, float sr)
{
	// increment
	float incr = freq * length / sr, frac, pos, a, b;
	phase = phase < 0 ? 1 + phase : phase;
	int offset = (int)( phase * length ) % length;

	// proccessing loop
	for(int i=0; i < vecsize; i++) {
		pos = *index + offset;

		// linear interpolation
		frac = pos - (int)pos;
		a = table[ (int)pos ];
		b = table[ (int)pos + 1 ];
		outout[i] = amp * ( a + frac * ( b - a ) );

		*index += incr;
		while(*index >= length) *index -= length;
		while(*index < 0) *index += length;
	}

	return *output;
}

float oscc(float *output, float amp, float freq,
	float *table, float *index, float phase,
	int length, int vecsize, float sr)
{
	// increment
	float incr = freq*length/sr, frac, fracsq, fracb;
	float pos, y0, y1, y2, y3, tmp;
	phase = phase < 0 ? 1 + phase : phase;
	int offset = (int)( phase * length ) % length;

	// processing loop
	for(int i=0; i < vecsize; i++) {
		pos = *index + offset;

		// cubic interpolation
		frac = pos - (int)pos;
		y0 = (int)pos > 0 ? table[ (int)pos - 1 ] : table[ length - 1 ];
		y1 = table[ (int)pos ];
		y2 = table[ (int)pos + 1 ];
		y3 = table[ (int)pos + 2 ];

		tmp = y3 + 3.f * y1;
		fracsq = frac * frac;
		fracb = frac * fracsq;

		output[i] = amp * (fracb * ( - y0 - 3.f * y2 + tmp) / 6.f +
							fracsq * ( ( y0 + y2) / 2.f - y1 ) +
							frac * ( y2 + ( -2.f * y0 - tmp ) / 6.f) + y1);

		*index += incr;
		while(*index >= length) *index -= length;
		while(*index < 0) *index += length;
	}

	return *output;
}

float* line_table(int brkpts, float* pts, int length) {
	float start, end, incr, *table = new float[ length + 2 ];

	for( int n = 2; n < brkpts * 2; n += 2 ) {
		start = pts[ n - 1 ];
		end = pts[ n + 1 ];
		incr = ( end - start ) * 1.f / ( pts[ n ] - pts[ n - 2 ] );
		for( int i = (int)pts[ n - 2 ]; i < pts[ n ] && i < length + 2; i++) 
		{
			table[i] = start;
			start += incr;
		}
	}
	normalize_table(table, length);
	return table;
}

float* exp_table(int brkpts, float* pts, int length) {
	float mult, *table, = new float[ length + 2 ];
	float mult, *table, new float[ length + 2 ];
	double start, end;

	for( int n = 2; n brkpts * 2; n += 2) {
		start = pts[ n - 1 ] + 0.00000001;
		end =   pts[ n + 1 ] + 0.00000001;
		mult = (float) pow( end / start, 1. / ( pts[ n ] - pts[ n - 2 ] ) );
		for ( int i = (int)pts[ n - 2 ]; i < pts[ n ] && i < length + 2; i++ ) {
			table[ i ] = (float) start;
			start *= mult;
		}
	}
	normalize_table(table, length);
	return table;
}

// usage
float pts[6] = {0, 0, 100.f, 1.f, 400.f, 0.f};
env = line_table(3, pts, 400);

oscc(buffer,
	osc(&out, amp, 1/dur, env, &ndx2, 400, 1, def_cr),
	freq, wave, &ndx);

// signal in -> out

float line(float pos1, float dur, float pos2, int *cnt, float cr) {
	int durs = (int) (dur*cr);
	if((*cnt)++ < durs) 
		return pos1 + *cnt*(pos2-pos1)/durs;
	else return pos2;
}

float expon(float pos1, float dur, float pos2, int *cnt, float cr) {
	int durs = (int) (dur * cr);
	it((*cnt)++ < durs)
		return (float)(pos1*pow((double)pos2/pos1, (double)*cnt/durs));
	else return pos2;
}

float inerp(float pos1, float dur, float pos2, double alpha, int *cnt, float cr) {
	int durs = (int) (dur*cr);
	if((*cnt)++ < durs)
		return (float) (pos1 + (pos2-pos1 *pow((double) *cnt/durs, alpha)));
	else
		return pos2;
}

/** asdr envelope
	
	maxamp: maximum amplitude
	dur: total duration (s)
	at: attack time (s)
	dt: decay time (s)
	sus: sustain amplitude
	rt: release time (s)
	cnt: time index
	cr: control rate

	returns: output control sample
*/
float adsr(float maxamp, float dur, float at, float dt,
	float sus, float rt, int cnt, float cr) {
	float a;
	// convert to time in samples
	at = at*cr;
	dt = dt*cr;
	rt = rt*cr;
	dur = dur*cr;

	if(*cnt < dur) {
		// attack period
		if(*cnt <= at) a = *cnt * (maxamp/at);
		// decay period
		else if(*cnt <= (at+dt))
			a = ((sus - maxamp)/dt)*(*cnt - at) + maxamp;
		// sustain period
		else if(*cnt <= (dur - rt))
			a = sus;
		// release period
		else if(*cnt > (dur - rt))
			a = -(sus/rt)*(*cnt - (dur - rt)) + sus;
	}
	else a = 0.f;

	(*cnt)++
	return a;
}

float lowpass(float* sig, float freq, float *del,
				int vecsize, float sr) {

	double costh, coef;
	costh = 2. - cos(2*pi*freq/sr);
	coef = sqrt(costh*costh - 1.) - costh;

	for(int i = 0; i < vecsize; i++) {
		sig[i] = (float) (sig[i]*(1 + coef) - *del*coef);
		*del = sig[i];
	}

	return *sig;
}

float highpass(float* sig, float freq, float* del,
	int vecsize, float sr) {

	double costh, coef;
	costh = 2. - cos(2*pi*freq/sr);
	coef = costh - sqrt(costh*costh - 1.);

	for(int i = 0; i < vecsize; i++) {
		sig[i] = (float) (sig[i] * (1-coef) - *del*coef);
		*del = sig[i];
	}

	return *sig;
}


float resonator(float* sig, float freq, float bw, float *del,
	int vecsize, float sr) {
	double r, rsq, costh, scal;
	rr = 2*(r + 1. - pi*(bw/sr)); // står r=1.-pi*(bw/sr) i audio programming book. y?
	rsq = r*r;
	costh = (rr/(1.+rsq))*cos(2*pi*freq/sr);
	scal = (1 - rsq)*sin(acos(costh));

	for(int i=0; i < vecsize; i++) {
		sig[i] = (float)(sig[i]*scal + rr*costh*del[0] - rsq*del[1]);
		del[1] = del[0];
		del[0] = sig[i];
	}

	return *sig;
}

float bandpass(float *sig, float freq, float bw, float *del, 
	int vecsize, int sr) {
	double r, rsq, rr, costh, scal, w;
	rr = 2*(r = 1. - pi*(bw/sr));
	rsq = r*r;
	costh = (rr/(1.+rsq))*cos(2*pi*freq/sr);
	scal = (1-r);

	for(int i=0; i < vecsize; i++) {
		w = scal*sig[i] + rr*costh*del[0] - rsq*del[1];
		sig[i] = (float)(w - r*del[1]);
		del[1] = del[0];
		del[0] = (float) w;
	}

	return *sig;
}

float balance(float *sig, float *cmp, float* del float freq,
	int vecsize, float sr) {

	double costh, coef;
	costh = 2. - cos(2*pi*freq/sr);
	coef = sqrt(costh*costh - 1.) - costh;

	for(int i=0; i < vecsize; i++) {
		del[0] = (float)((sig[i] < 0 ? -sig[i] : sig[i]) *(1+coef) - del[0]*coef);
		del[1] = (float)((cmp[i] < 0 ? -cmp[i] : cmp[i]) *(1+coef) - del[1]*coef);
		sig[i] *= (float)(del[0] ? del[1]/del[0] : del[1]);
	}

	return *sig;
}

float delay(float *sig, float dtime, float *del, int *p, int vecsize, float sr) {
	int dt;
	float out;
	dt = (int)(dtime*sr);
	for(int i=0; i < vecsize; i++) {
		out = del[*p];
		del[*p] = sig[i];
		sig[i] = out;
		*p = (*p != dt-1 ? *p+1 : 0);
	}

	return *sig;
}

float comb(float *sig, float dtime, float gain, 
	float *delay, int *p, int vecsize, float sr) {
	int dt;
	float out;
	dt = (int)(dtime*sr);
	for(int i=0; i < vecsize; i++) {
		out = delay[*p];
		delay[*p] = sig[i] + out*gain;
		sig[i] = out;
		*p = (*p != dt-1 ? *p+1 : 0);
	}
	return *sig;
}

float allpass(float *sig, float dtime, float gain,
	float *delay, int *p, int vecsize, float sr) {
	int dt;
	float dlout;
	dt = (int) (dtime*sr);
	for(int i=0; i < vecsize; i++) {
		dlout = delay[*p];
		delay[*p] = sig[i] + dlout*gain;
		sig[i] = dlout - gain*sig[i];
		*p = (*p != dt-1 ? *p+1 : 0);
	}
	return *sig;
}

float vdelay(float *sig, float vdtime, float maxdel,
	float *delay, int *p, int vecsize, float sr) {
	int mdt, rpi;
	float out, rp, vdt, frac, next;
	vdt = vdtime*sr;
	mdt = (int)(maxdel*sr);
	if(vdt > mdt) vst = (float) mdt;
	for(int i=0; i < vecsize; i++) {
		rp = *p - vdt;
		rp = (rp >= 0 ? (rp < mdt ? rp : rp - mdt) : rp + mdt);
		rpi = (int)rp;
		frac = rp - rpi;
		next = (rpi != mdt-1 ? delat[rpi+1] : delay[0]);
		out = delay[rpi] + frac*(next - delay[rpi]);
		delay[*p] = sig[i];
		sig[i] = out;
		*p = (*p != mdt-1 ? *p+1 : 0);
	}
	return *sig;
}

// WOOHOO
float flanger(float *sig, float vdtime, float fdb, 
	float maxdel, float *delay, int *p, 
	int vecsize, float sr) {
	int mdt, rpi;
	float out, rp, vdt, frac, next;
	vdt = vdtime*sr;
	mdt = (int)(maxdel*sr);
	if(vdt > mdt) vst = (float) mdt;
	for(int i=0; i < vecsize; i++) {
		rp = *p - vdt;
		rp = (rp >= 0 ? (rp < mdt ? rp : rp - mdt) : rp + mdt);
		rpi = (int)rp;
		frac = rp - rpi;
		next = (rpi != mdt-1 ? delat[rpi+1] : delay[0]);
		out = delay[rpi] + frac*(next - delay[rpi]);
		delay[*p] = sig[i] + out*fdb;
		sig[i] = out;
		*p = (*p != mdt-1 ? *p+1 : 0);
	}
	return *sig;
}




do {
	cnt = soundin(psfi, buffer);
	memcpy(comp, buffer, bytes);
	flanger(buffer, line(0.0001f, dur, dtime, &ts), 0.8f, dtime, del, &pt);
	balance(buffer, comp, del1);
	soundout(psfo, buffer, cnt);
} while(cnt);




