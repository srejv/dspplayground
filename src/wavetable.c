
#include "wavetable.h"

WaveTable* WaveTable_create() {
	WaveTable* table = (WaveTable*)malloc(sizeof(WaveTable));

	// non inverting unless user sets it
	table->bInvert = false;

	// slope and y-intercept values for the
	// Triangle Wave
	// rising edge1:
	float mt1 = 1.0/256.0;
	float bt1 = 0.0;

	// rising edge 2:
	float mt2 = 1.0/256.0;
	float bt2 = -1.0;

	// valling edge:
	float mtf2 = -2.0/512.0;
	float btf2 = 1.0;

	// Sawtooth
	// rising edge1:
	float ms1 = 1.0/512.0;
	float bs1 = 0.0;

	// rising edge2:
	float ms2 = 1.0/512.0;
	float bs2 = -1.0;

	float fMaxTri = 0;
	float fMaxSaw = 0;
	float fMaxSqr = 0;

	// setup arrays
	for(int i = 0; i < 1024; i++) 
	{
		// sample the sinusoid, 1024 points
		// sin(wnT) = sin(2pi*i/1024)
		table->SinArray[i] = sin( ((float)i/1024.0) * (2*pi) );

		// saw, triangle and square are just logic mechanisms
		// Sawtooth
		table->SawtoothArray[i] = i < 512 ? ms1*i + bs1 : ms2*(i-511) + bs2;

		// Triangle
		if( i < 256 )
			table->TriangleArray[i] = mt1*i + bt1; // mx + b; rising edge 1
		else if(i < 768)
			table->TriangleArray[i] = mtf2*(i-256) + btf2; // mx + b; falling edge
		else 
			table->TriangleArray[i] = mt2*(i-768) + bt2; // mx + b; rising edge 2

		// square
		if(i!=0)
			table->SquareArray[i] = i < 512 ? +1.0 : -1.0;
		else 
			table->SquareArray[i] = 0;

		// zero to start, then loops build the rest
		table->SawtoothArray_BL5[i] = 0.0;
		table->SquareArray_BL5[i] = 0.0;
		table->TriangleArray_BL5[i] = 0.0;

		// sawtooth: += (-1)^g+1 * (1/g) * sin(wnT)
		for(int g=1; g <= 6; g++)
		{
			double n = double(g);
			table->SawtoothArray_BL5[i] += pow((float)-1.0, (float)(g+1))*(1.0/n)*sin(2.0*pi*i*n/1024.0);
			// down saw table->SawtoothArray_BL5[i] += (1.0/n) * sin(2.0*pi*i*n/1024.0);
		}

		// triangle: += (-1)^g(1/(2g+1+^2))sin(w(2n+1)T)
		for(int g=0; g<=3; g++)
		{
			double n = double(g);
			table->TriangleArray_BL5[i] += pow((float)-1.0, (float)n)*(1.0/pow((float)(2*n + 1), (float)2.0))*sin(2.0*pi*(2.0*n + 1)*i/1024.0);
		}

		// square: += (1/g)sin(wnT)
		for(int g=1; g<8; g+=2)
		{
			double n = double(g);
			table->SquareArray_BL5[i] += (1.0/n)*sin(2.0*pi*i*n/1024.0);
		}

		// store the max values
		if(i == 0)
		{
			fMaxSaw = m_SawtoothArray_BL5[i];
			fMaxTri = m_TriangleArray_BL5[i];
			fMaxSqr = m_SquareArray_BL5[i];
		}
		else
		{
			// test and store
			if(table->SawtoothArray_BL5[i] > fMaxSaw)
				fMaxSaw = table->m_SawtoothArray_BL5[i];

			if(table->TriangleArray_BL5[i] > fMaxTri)
				fMaxTri = table->m_TriangleArray_BL5[i];

			if(table->SquareArray_BL5[i] > fMaxSqr)
				fMaxSqr = table->SquareArray_BL5[i];
		}
	}

	// normalize the bandlimited tables
	for(int i = 0; i < 1024; i++)
	{
		// normalize it
		table->SawtoothArray_BL5[i] /= fMaxSaw;
		table->TriangleArray_BL5[i] /= fMaxTri;
		table->SquareArray_BL5[i] /= fMaxSqr;
	}

	// clear variables
	table->fReadIndex = 0.0;
	table->fQuadPhaseReadIndex = 0.0;
	table->fInc = 0.0;

	// initialize inc
	WaveTable_reset(table);
	table->fFrequency_Hz = 440;
	table->uOscType = 0;
	table->uTableMode = 0;
	table->uPolarity = 0;

	WaveTable_cookFrequency(table);

	return table;
}


void WaveTable_prepareForPlay(WaveTable *this) {
	// reset the index
	WaveTable_reset(this);

	// cook curent frequency
	WaveTable_cookFrequency(this);
}
void WaveTable_doOscillate(WaveTable* this, float* pYn, float* pYqn)
{
	// our output value for this cycle
	float fOutSample = 0;
	float fQuadPhaseOutSample = 0;

	// get INT part
	int nReadIndex = (int)this->fReadIndex;
	int nQuadPhaseReadIndex = (int)this->fQuadPhaseReadIndex;

	// get FRAC part
	float fFrac = this->fReadIndex - nReadIndex;
	float fQuadFrac = this->fQuadPhaseReadIndex - nQuadPhaseReadIndex;

	// setup second index for interpolation; wrap the buffer if needed
	int nReadIndexNext = nReadIndex + 1 > 1023 ? 0 :  nReadIndex + 1;
	int nQuadPhaseReadIndexNext = nQuadPhaseReadIndex + 1 > 1023 ? 0 :  nQuadPhaseReadIndex + 1;

	// interpolate the output
	switch(this->uOscType)
	{
		case sine:
			fOutSample = dLinTerp(0, 1, this->SinArray[nReadIndex], this->SinArray[nReadIndexNext], fFrac);
			fQuadPhaseOutSample = dLinTerp(0, 1, this->SinArray[nQuadPhaseReadIndex], this->SinArray[nQuadPhaseReadIndexNext], fQuadFrac);
			break;

		case saw:
			if(this->uTableMode == normal)	// normal
			{
				fOutSample = dLinTerp(0, 1, this->SawtoothArray[nReadIndex], m_SawtoothArray[nReadIndexNext], fFrac);
				fQuadPhaseOutSample = dLinTerp(0, 1, this->SawtoothArray[nQuadPhaseReadIndex], m_SawtoothArray[nQuadPhaseReadIndexNext], fQuadFrac);
			}
			else						// bandlimited
			{
				fOutSample = dLinTerp(0, 1, this->SawtoothArray_BL5[nReadIndex], this->SawtoothArray_BL5[nReadIndexNext], fFrac);
				fQuadPhaseOutSample = dLinTerp(0, 1, this->SawtoothArray_BL5[nQuadPhaseReadIndex], this->SawtoothArray_BL5[nQuadPhaseReadIndexNext], fQuadFrac);
			}
			break;

		case tri:
			if(this->uTableMode == normal)	// normal
			{
				fOutSample = dLinTerp(0, 1, this->TriangleArray[nReadIndex], this->TriangleArray[nReadIndexNext], fFrac);
				fQuadPhaseOutSample = dLinTerp(0, 1, this->TriangleArray[nQuadPhaseReadIndex], this->TriangleArray[nQuadPhaseReadIndexNext], fQuadFrac);
			}
			else						// bandlimited
			{
				fOutSample = dLinTerp(0, 1, this->TriangleArray_BL5[nReadIndex], this->TriangleArray_BL5[nReadIndexNext], fFrac);
				fQuadPhaseOutSample = dLinTerp(0, 1, this->TriangleArray_BL5[nQuadPhaseReadIndex], this->TriangleArray_BL5[nQuadPhaseReadIndexNext], fQuadFrac);
			}
			break;

		case square:
			if(this->uTableMode == normal)	// normal
			{
				fOutSample = dLinTerp(0, 1, this->SquareArray[nReadIndex], this->SquareArray[nReadIndexNext], fFrac);
				fQuadPhaseOutSample = dLinTerp(0, 1, this->SquareArray[nQuadPhaseReadIndex], this->SquareArray[nQuadPhaseReadIndexNext], fQuadFrac);
			}
			else						// bandlimited
			{
				fOutSample = dLinTerp(0, 1, this->SquareArray_BL5[nReadIndex], this->SquareArray_BL5[nReadIndexNext], fFrac);
				fQuadPhaseOutSample = dLinTerp(0, 1, this->SquareArray_BL5[nQuadPhaseReadIndex], this->SquareArray_BL5[nQuadPhaseReadIndexNext], fQuadFrac);
			}
			break;

		// always need a default
		default:
			fOutSample = dLinTerp(0, 1, this->SinArray[nReadIndex], this->SinArray[nReadIndexNext], fFrac);
			fQuadPhaseOutSample = dLinTerp(0, 1, this->SinArray[nQuadPhaseReadIndex], this->SinArray[nQuadPhaseReadIndexNext], fQuadFrac);
			break;
	}

	// add the increment for next time
	this->fReadIndex += this->fInc;
	this->fQuadPhaseReadIndex += this->fInc;

	// check for wrap
	if(this->fReadIndex >= 1024)
		this->fReadIndex = this->fReadIndex - 1024;

	if(this->fQuadPhaseReadIndex >= 1024)
		this->fQuadPhaseReadIndex = this->fQuadPhaseReadIndex - 1024;

	// write out
	*pYn = fOutSample;
	*pYqn = fQuadPhaseOutSample;

	if(m_bInvert)
	{
		*pYn *= -1.0;
		*pYqn *= -1.0;
	}

	if(m_uPolarity == unipolar)
	{
		*pYn *= 0.5;
		*pYn += 0.5;

		*pYqn *= 0.5;
		*pYqn += 0.5;
	}
}

void WaveTable_reset(WaveTable* this) {
	table->fReadIndex = 0;
	table->fQuadPhaseReadIndex = 0;
}

void WaveTable_setSampleRate(WaveTable* this, float f) {
	this->nSampleRate = nSampleRate;
	WaveTable_cookFrequency(this);
}

void WaveTable_cookFrequency(WaveTable* this) {
	// inc = L*fd/fs
	this->fInc = 1024.0*this->fFrequency_Hz/(float)this->nSampleRate;
}
