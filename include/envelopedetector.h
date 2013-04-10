/**
	EnvelopeDetector.

	Modes: PEAK = 0, MS = 1, RMS = 2

	Modified version from musicdsp.org (originally a C++ class)
*/

#ifndef ENVELOPEDETECTOR_H_
#define ENVELOPEDETECTOR_H_

// har att göra med övergången (för digitala signaler är det vid 1%, analoga 36.7% eller något sånt)
#define DIGITAL_TC -2.0; // log(1%)
#define ANALOG_TC -0.43533393574791066201247090699309; // log(36.7%)

typedef struct EnvelopeDetector_t {
	int nSample;
	float fAttackTime;
	float fReleaseTime;
	float fAttackTime_mSec;
	float fReleaseTime_mSec;
	float fSampleRate;
	float fEnvelope;
	unsigned int uDetectMode;
	bool bAnalogTC;
	bool bLogDetector;
} EnvelopeDetector;

EnvelopeDetector* EnvelopeDetector_create();
void EnvelopeDetector_init(EnvelopeDetector* detector, float samplerate, float attackInMs, float releaseInMs, bool bAnalogTC, unsigned int uDetect, bool bLogDetector);
// after init
void EnvelopeDetector_setTCModeAnalog(EnvelopeDetector* detector, bool bAnalogTC);
void EnvelopeDetector_setAttackTime(EnvelopeDetector* detector, float attackInMs);
void EnvelopeDetector_setReleaseTime(EnvelopeDetector* detector, float releaseInMs);

// Codes:
// DETECT PEAK 	= 0
// DETECT MS	= 1
// DETECT RMS	= 2
//
void EnvelopeDetector_setDetectMode(EnvelopeDetector* detector, unsigned int uDetect);

void EnvelopeDetector_setSampleRate(EnvelopeDetector* detector, float f);
void EnvelopeDetector_setLogDetect(EnvelopeDetector* detector, bool b);

// call this to detect; it returns the peak ms or rms value at that instant
float EnvelopeDetector_detect(EnvelopeDetector* detector, float fInput);

// call to reset the detector	
void EnvelopeDetector_reset(EnvelopeDetector* detector);

#endif