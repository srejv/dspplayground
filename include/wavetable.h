
#ifndef WAVETABLE_H_
#define WAVETABLE_H_

typedef struct WaveTable_t {
	float SinArray[1024];
	float SawtoothArray[1024];
	float TriangleArray[1024];
	float SquareArray[1024];

	// band limited to 5 partials (because of NyQuist)
	float SawtoothArray_BL5[1024];
	float TriangleArray_BL5[1024];
	float SquareArray_BL5[1024];

	// current read location
	float fReadIndex;
	float fQuadPhaseReadIndex;

	// our inc value
	float fInc;

	// fs value
	int nSampleRate;

	// user controlled variables
	// frequency
	float fFrequency_Hz;

	// inverted output
	bool bInvert;

	// type
	unsigned int uOscType;
	enum {sine, saw, tri, square};

	// mode
	unsigned int uTableMode;
	enum {normal, bandlimited};

	// polarity
	unsigned int uPolarity;
	enum { bipolar, unipolar };

} WaveTable;

WaveTable* WaveTable_create();

void WaveTable_prepareForPlay();
void WaveTable_doOscillate();

void WaveTable_reset();

void WaveTable_setSampleRate();
void WaveTable_cookFrequency();

#endif