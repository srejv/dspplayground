


#define Q 15
#define K ( 1 << ( Q - 1 ) )

// add
signed int a, b, result;
result = a + b;

// sub
ssigned int a, b, result;
result = a - b;

// mult
signed int a, b, result;
signed long int temp;
temp = (long int) a * (long int)  b; // result type is operand's type
// Rounding; mid values are rounded up
temp += K;
// Correct by dividing by base
result = temp >> Q;

// div
signed int a, b, result;
signed long int temp;
// pre-multiply by the base (Upscale to Q16 so that the result will be in Q8 format)
temp = (long int)a << Q;
// So the result will be rounded ; mid values are rounded up.
temp = temp + b/2;
result = temp / b

//
a/b = a * (1/b);




#define short s16
#define long s32

s16 x,y,z;
z = (s16)( ( (s32)x * (s32)y )  >> 15 );

// Använd inte -32768, blir oveflow.
[-32767, 32767]

// Addition produces one extra bit

#include <stdint.h>

static const int16_t Q15_N = 15; /* fractional bits, n in Q15, Qm.15 format */
#define Q15 int16_t
#define Q30	int32_t

Q15 Q15_add(Q15 a, Q15 b) {
	return a + b;
}

Q15 Q15_subtract(Q15 a, Q15 b) {
	return a - b;
}

Q15 Q15_multiply(Q15 a, Q15 b) {
	/* Rounding: mid values are rounded up */
	static const short Q15_ONE_HALF = 1 << (Q15_N - 1);
	Q30 resultTimes2N = (Q30)a * b + Q15_ONE_HALF;
	/* Correct by dividing by base */
	return (Q15)(resultTimes2N >> Q15_N);
}

Q15 Q15_divide(Q15 a, Q15 b) {
	/* pre-multiply by the base */
	Q30 aTimes2N = (Q30)a << Q15_N;
	/* So the result will be rounded; mid values are rounded up */
	Q30 roundedATimes2N = aTimes2N + b/2;
	return (Q15)(roundedATimes2N / b);
}