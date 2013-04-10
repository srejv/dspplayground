#include <stdlib.h>
#include <limits.h>

static double inv_randmax = 1.0 / RAND_MAX;

double trirand(void)
{
	double r1, r2;
	r1 = rand() * inv_randmax;
	r2 = rand() * inv_randmax;
	return (r1 + r2) * 0.5;
}

// usage
double tpdf;
float fsamp;
short samp;
tpdf = ( 4.0 * trirand() ) - 2.0;
samp = (short) ( fsamp * 32765.0 + tpdf );