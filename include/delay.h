
#ifndef ___INCLUDED_DELAY_H___ 
#define ___INCLUDED_DELAY_H___ 

typedef struct delay_t {
	double sr;
	double delay;
	int wp;	
	float gain;
	float feedback;
	int length;
	float samples[];
} DELAY;

/*
typedef struct table_t {
	unsigned int length;
	float samples[];
} TABLE;

TABLE * new_table() {
	// allocate memory
}

TABLE * new_sine_table() {

}
*/


/*
typedef struct vdelay_t {
	double sr;
	double delay;
	int wp;	  // write pointer
	float rp; // read pointer

	float gain;
	float feedback;

	float input();
	float output();
} VDELAY;
*/

DELAY * new_delay(float delay_time, float sample_rate, float gain, float feedback);
void delay_free(DELAY** d);
void delay_tick(DELAY* d);
float delay_processframe(DELAY * self, float in);
float delay_processblock(DELAY * self, float *in, float *out, unsigned long nFrames, int nChannels);

/*
VDELAY * new_vdelay(float delay_time, float sample_rate, float gain, float feedback);
void free_vdelay(VDELAY** self);
void vdelay_tick(VDELAY* self)
float vdelay_process_sample(VDELAY * self, float in);
*/
#endif