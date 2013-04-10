RDelay* RDelay_create(int nDelayLength)
{

	RDelay* r = (RDelay*)malloc(sizeof(RDelay) + sizeof(float)*this->nBufferSize)
	r->pBuffer = NULL;

	r->fOutputAttenuation_dB = 0;
	r->fDelay_ms = 0.0;

	r->fOutputAttenuation = 0.0;
	r->fDelayInSamples = 0.0;
	r->nSampleRate = 0;

	resetDelay(r);

	return r;
}

void RDelay_destroy(RDelay* this)
{
	if(this) {
		if(this->pBuffer)
			free(this->pBuffer);

		this->pBuffer = NULL;

		free(this);
		this = NULL;
	}
}

void RDelay_init(RDelay* this, int nDelayLength)
{
	this->nBufferSize = nDelayLength;
	//this->pBuffer = (float*)malloc(sizeof(float)*this->nBufferSize]);
	// flush buffer
	memset(this->pBuffer, 0, this->nBufferSize*sizeof(float));
}

void RDelay_resetDelay(RDelay* this)
{
	// flush buffer
	if(this->pBuffer)
		memset(this->pBuffer, 0, this->nBufferSize*sizeof(float));

	// init read/write indices
	this->nWriteIndex = 0; // reset the Write index to top
	this->nReadIndex = 0; // reset the Write index to top

	RDelay_cookVariables(this);
}

void RDelay_setDelay_mSec(RDelay* this, float fmSec)
{
	this->fDelay_ms = fmSec;
	RDelay_cookVariables(this);
}

void RDelay_setOutputAttenuation_dB(RDelay* this, float fAttendB)
{
	this->fOutputAttenuation_dB = fAttendB;
	RDelay_cookVariables(this);
}

void RDelay_cookVariables(RDelay* this)
{
	this->fOutputAttenuation = pow((float)10.0, (float)this->fOutputAttenuation_dB/(float)20.0);

	this->fDelayInSamples = this->fDelay_ms*((float)this->nSampleRate/1000.0);

	// subtract to make read index
	this->nReadIndex = this->nWriteIndex - (int)this->fDelayInSamples;

	//  the check and wrap BACKWARDS if the index is negative
	if (this->nReadIndex < 0)
		this->nReadIndex += this->nBufferSize;	// amount of wrap is Read + Length

}

void RDelay_writeDelayAndInc(RDelay* this, float fDelayInput)
{
	// write to the delay line
	this->pBuffer[this->nWriteIndex] = fDelayInput; // external feedback sample

	// incremnent the pointers and wrap if necessary
	this->nWriteIndex++;
	if(this->nWriteIndex >= this->nBufferSize)
		this->nWriteIndex = 0;

	this->nReadIndex++;
	if(this->nReadIndex >= this->nBufferSize)
		this->nReadIndex = 0;
}

float RDelay_readDelay(RDelay* this)
{
	// Read the output of the delay at m_nReadIndex
	float yn = this->pBuffer[this->nReadIndex];

	// Read the location ONE BEHIND yn at y(n-1)
	int nReadIndex_1 = this->nReadIndex - 1;
	if(nReadIndex_1 < 0)
		nReadIndex_1 = this->nBufferSize-1; // m_nBufferSize-1 is last location

	// get y(n-1)
	float yn_1 = this->pBuffer[nReadIndex_1];

	// interpolate: (0, yn) and (1, yn_1) by the amount fracDelay
	float fFracDelay = this->fDelayInSamples - (int)this->fDelayInSamples;

	return dLinTerp(0, 1, yn, yn_1, fFracDelay); // interp frac between them
}

float RDelay_readDelayAt(RDelay* this, float fmSec)
{
	float fDelayInSamples = fmSec*((float)this->nSampleRate)/1000.0;

	// subtract to make read index
	int nReadIndex = this->nWriteIndex - (int)fDelayInSamples;

	// Read the output of the delay at m_nReadIndex
	float yn = this->pBuffer[nReadIndex];

	// Read the location ONE BEHIND yn at y(n-1)
	int nReadIndex_1 = nReadIndex - 1;
	if(nReadIndex_1 < 0)
		nReadIndex_1 = this->nBufferSize-1; // m_nBufferSize-1 is last location

	// get y(n-1)
	float yn_1 = this->pBuffer[nReadIndex_1];

	// interpolate: (0, yn) and (1, yn_1) by the amount fracDelay
	float fFracDelay = fDelayInSamples - (int)fDelayInSamples;

	return dLinTerp(0, 1, yn, yn_1, fFracDelay); // interp frac between them
}

bool RDelay_processAudio(RDelay* this, float* pInput, float* pOutput)
{
	// Read the Input
	float xn = *pInput;

	// read delayed output
	float yn = this->fDelayInSamples == 0 ? xn : RDelay_readDelay(this);

	// write to the delay line
	RDelay_writeDelayAndInc(this, xn);

	// output attenuation
	*pOutput = this->fOutputAttenuation*yn;

	return true; // all OK
}
