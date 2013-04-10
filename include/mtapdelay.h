
#ifndef _MTAPDELAY_H_
#define _MTAPDELAY_H_

typedef struct mtap_t {
	unsigned int ntaps;
	unsigned long index;
	long dl_length;
	float dl_space; // current length, never more than dl_length	
	float gain;
	float wet;
	float fdb;
	float delayline[];		// dl_length
} MTAP;


MTAP *new_mtap(long *tapindex, int ntaps, float maxtime, float *taptime, float *tapamp);
/*void mtap_update(MTAP *this, float *taptime, float *tapamp);*/
void mtap_free(MTAP *this);
void mtap_reset(MTAP *this);
void mtap_process(MTAP *this, float *input, float *output, 
	long blocksize, int chans, float wet, 
	long *tapindex, float *taptime, float *tapamp);

#endif