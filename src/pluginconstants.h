
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// this #define enables teh following constants from math.h
#define _MATH_DEFINES_DEFINED
/* 
#define M_E			2.718...
#define M_LOG2E		1,442695...
#define M_LOG10		0.4342944819...	
#define M_LN2		9,68
#define M_LN10		
#define M_PI		
#define M_PI_2		
#define M_PI_4		
#define M_1_PI		
#define M_2_SQRTPI	
#define M_SQRT2		
#define m_SQRT1_2	
*/
#include <math.h>

#define FLT_EPSILON_PLUS	 1.192092896e-07		// smallest such that 1.0+FLT_EPSILON != 1.0 
#define FLT_EPSILON_MINUS	-1.192092896e-07		// smallest such that 1.0-FLT_EPSILON != 1.0
#define FLT_MIN_PLUS		 1.175494351e-38		// min positive value
#define FLT_MIN_MINUS		-1.175494351e-38		// min negative values

/*
	Function:	lagrpol() implements n-order Lagrange Interpolation

	Inputs:		double* x	Pointer to an array containing x-coordinates of the input values
				double* y	Pointer to an array containing y-coordinates of the input values
				int n		The order of the interpolator, this is also the length of the x,y input arrays
				double xbar	The x-coorinates whose y-value we want to inerpolate

	Returns		The interpolated value y at xbar. xbar ideally is between the middle two values in the input array,
				but can be anywhere within the limits, which is needed for interpolating the first few or last few samples
				in a table with a fixed size
*/
inline double lagrpol(double* x, double* y, int n, double xbar)
{
	int i,j;
	double fx=0.0;
	double l =1.0;
	for(i=0; i<n; i++)
	{
		l=1.0;
		for(j=0; j<n; j++)
		{
			if(j!=i)
				l *= (xbar-x[j])/(x[i]-x[j]);
		}
		fx += l*y[i];
	}
	return fx;
}

inline float dLinTerp(float x1, float x2, float y1, float y2, float x)
{
	float denom = x2 - x1;
	if(denom == 0)
		return y1; // should never happen

	// calculate decimal position of x;
	float dx = (x - x1)/(x2 - x1);

	// use weighted sum method of interpolating
	float result = dx*y2 + (1-dx)*y1;

	return result;
}

inline bool normalizeBuffer(double* pInputBuffer, unsigned int uBufferSize)
{
	double fMax = 0;

	for(unsigned int j=0; j < uBufferSize; j++)
	{
		if((fabs(pInputBuffer[j])) > fMax)
			fMax = fabs(pInputBuffer[j]);
	}


	if(fMax > 0)
	{
		fMax = 1.0/fMax;
		for(unsigned int j=0; j < uBufferSize; j++)
			pInputBuffer[j] = pInputBuffer[j]*fMax;
	}

	return true;
}

