
#include "envelopedetector.h"


EnvelopeDetector* EnvelopeDetector_create() {
	EnvelopeDetector* detector = (EnvelopeDetector*)malloc(sizeof(EnvelopeDetector));
	detector->fAttackTime_mSec = 0.0;
	detector->fReleaseTime_mSec = 0.0;
	detector->fAttackTime = 0.0;
	detector->fReleaseTime = 0.0;
	detector->fSampleRate = 44100;
	detector->fEnvelope = 0.0;
	detector->uDetectMode = 0;
	detector->nSample = 0;
	detector->bAnalogTC = false;
	detector->bLogDetector = false;
}

void EnvelopeDetector_reset(EnvelopeDetector* detector) {
	detector->fEnvelope = 0.0;
	detector->nSample = 0;
}

void EnvelopeDetector_init(EnvelopeDetector* detector, float samplerate, float attackInMs, float releaseInMs, bool bAnalogTC, unsigned int uDetect, bool bLogDetector)
{
	detector->fEnvelope = 0.0;
	detector->fSampleRate = samplerate;
	detector->bAnalogTC = bAnalogTC;
	detector->fAttackTime_mSec = attackInMs;
	detector->fReleaseTime_mSec = releaseInMs;
	detector->uDetectMode = uDetect;
	detector->bLogDetector = bLogDetector;

	// set them
	EnvelopeDetector_setAttackTime(detector, attackInMs);
	EnvelopeDetector_setReleaseTime(detector, releaseInMs);
}

void EnvelopeDetector_setAttackTime(EnvelopeDetector* detector, float attackInMs)
{
	detector->fAttackTime_mSec = attackInMs;

	if(detector->bAnalogTC)
		detector->fAttackTime = exp(ANALOG_TC / ( attackInMs * detector->fSampleRate * 0.001 ) );
	else
		detector->fAttackTime = exp(DIGITAL_TC / ( attackInMs * detector->fSampleRate * 0.001 ) );
}

void EnvelopeDetector_setReleaseTime(EnvelopeDetector* detector, float releaseInMs)
{
	detector->fReleaseTime_mSec = releaseInMs;

	if(detector->bAnalogTC)
		detector->fReleaseTime = exp(ANALOG_TC / ( releaseInMs * detector->fSampleRate * 0.001 ) );
	else
		detector->fReleaseTime = exp(DIGITAL_TC / ( releaseInMs * detector->fSampleRate * 0.001 ) );
}

void EnvelopeDetector_setTCModeAnalog(EnvelopeDetector* detector, bool bAnalogTC)
{
	detector->bAnalogTC = bAnalogTC;
	EnvelopeDetector_setAttackTime(detector, detector->fAttackTime_mSec);
	EnvelopeDetector_setReleaseTime(detector, detector->fReleaseTime_mSec);
}

void EnvelopeDetector_setDetectMode(EnvelopeDetector* this, unsigned int uDetect) {
	this->uDetectMode = uDetectMode;
}

void EnvelopeDetector_setSampleRate(EnvelopeDetector* this, float f) {
	this->fSampleRate = f;
}

void EnvelopeDetector_setLogDetect(EnvelopeDetector* this, bool b) {
	this->bLogDetector = b;
}

void EnvelopeDetector_detect(EnvelopeDetector* this, float fInput)
{
	switch(this->uDetectMode)
	{
	case 0:
		fInput = fabs(fInput);
		break;
	case 1:
		fInput = fabs(fInput) * fabs(input);
		break;
	case 2:
		fInput = pow((float)fabs(fInput) * (float)fabs(fInput), (float)0.5);
		break;
	default:
		fInput = fabs(fInput);
		break;
	}

	float fOld = this->fEnvelope;
	if(fInput > this->fEnvelope)
		this->fEnvelope = this->fAttackTime * (this->fEnvelope - fInput) + fInput;
	else
		this->fEnvelope = this->fReleaseTime * (this->fEnvelope - fInput) + fInput;

	if(this->fEnvelope > 0.0 && this->fEnvelope < FLT_MIN_PLUS) this->fEnvelope = 0;
	if(this->fEnvelope < 0.0 && this->fEnvelope > FLT_MIN_MINUS) this->fEnvelope = 0;

	// bound them; can happen when using pre-detector gains of more than 1.0
	this->fEnvelope = min(this->fEnvelope, 1.0);
	this->fEnvelope = max(this->fEnvelope, 0.0);

	// 16-bit scaling!
	if(this->bLogDetector)
	{
		if(this->fEnvelope <= 0)
			return -96.0; // 16 bit noise floor

		return 20*log10(this->fEnvelope);
	}

	return this->fEnvelope;
}


