#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "euler_utils.h"

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define LOG_GR_REC 2.0780869212350273

/**
 ********************************************************
 * Various utility functions for mathematical problems. 
 * Author: Hjalti Thor Isleifsson
 ********************************************************
 */
int gcd(int a, int b) {
	int r,q;
	while (b != 0) {
		q = a / b;
		r = a - q * b;
		a = b;
		b = r;
	}
	return a;
}

int lcm(int a, int b) {
	return (a * b) / gcd(a, b);
}


void bezout(int m, int n, int* a, int* b, int* c) {
	int q = m / n;
	int old_r = n;
	int r = m - q * n;

	int old_a = 1;
	int ma = 0;
	int old_b = 0;
	int mb = 1;

	int old_old_a, old_old_b, old_old_r;

	while (r != 0) {
		old_old_a = old_a;
		old_a = ma;
		ma = old_old_a - q * ma;

		old_old_b = old_b;
		old_b = mb;
		mb = old_old_b - q * mb;

		q = old_r / r;
		
		old_old_r = old_r;
		old_r = r;
		r = old_old_r - q * r;
	}

	*a = ma;
	*b = mb;
	*c = old_r;
}

int inv(int x, int p) {
	int a, b, c;
	bezout(x, p, &a, &b, &c);
	if (c != 1) {
		printf("%d does not have an inverse in Z_%d", x,p);
		return 0;
	} else if (a < 0) {
		return a + p;
	} else {
	 	return a;
	}
}

bool isprime(int n) {
	if (n < 4) {
		return n > 1;
	} else if ((n & 1) == 0 || (n % 3) == 0) {
		return false;
	} else {
		int d = 5;
		while (d * d <= n) {
			if (n % d == 0 || n % (d + 2) == 0) {
				return false;
			}
			d += 6;
		}
		return true;
	}
}

//All primes <= n
void allprimes(int n, int** primes, int* n_primes) {
	if (n < 2) {
		int* mprimes = (int*)malloc(0);
		*primes = mprimes;
		*n_primes = 0;
		return;
	}
	if (n < 3) {
		int* mprimes = (int*)malloc(1 * sizeof(int));
		mprimes[0] = 2;
		*primes = mprimes;
		*n_primes = 1;
		return;
	}
	if (n < 5) {
		int* mprimes = (int*)malloc(2 * sizeof(int));
		mprimes[0] = 2;
		mprimes[1] = 3;
		*primes = mprimes;
		*n_primes = 2;
		return;
	}

	const int mp_len = round(1.2 * n / log(n));
	const int d = isqrt(n); //Interval length.

	int* mprimes = (int *)malloc(mp_len * sizeof(int));
	bool flags[d];

	int pc = 0; //Prime counter
	int mp_idx, f_idx, f_idx2;
	
	for (mp_idx = 0; mp_idx < mp_len; ++mp_idx) {
		mprimes[mp_idx] = 0;
	}

	flags[0] = false;
	flags[1] = false;
	for (f_idx = 2; f_idx < d; ++f_idx) { 
		flags[f_idx] = true;
	}

	for (f_idx = 2; f_idx < d; ++f_idx) {
		if (flags[f_idx]) {
			mprimes[pc++] = f_idx;
			for (f_idx2 = 2 * f_idx; f_idx2 < d; f_idx2 += f_idx) {
				flags[f_idx2] = false;
			}
		}
	}

	const int n_iv = (n + d - 1) / d; //Number of intervals.

	int offset, upto, sqrt_upto, p, f_upto, iv_idx;
	for (iv_idx = 1; iv_idx < n_iv; ++iv_idx) {
		//Must clean flags
		for (f_idx = 0; f_idx < d; ++f_idx) { flags[f_idx] = true; }

		offset = iv_idx * d;
		upto = MIN(n + 1, offset + d);

		sqrt_upto = isqrt(upto);
		f_upto = upto - offset;

		for (mp_idx = 0; mp_idx < mp_len; ++mp_idx) {
			p = mprimes[mp_idx];
			if (p > sqrt_upto || p == 0) {
				break;
			} else {
				for (f_idx = ((offset + p - 1) / p) * p - offset; f_idx < f_upto; f_idx += p) {
					flags[f_idx] = false;
				}
			}
		}

		for (f_idx = 0; f_idx < f_upto; ++f_idx) {
			if (flags[f_idx]) {
				mprimes[pc++] = f_idx + offset;
			}
		}
	}
	*n_primes = pc;
	mprimes = (int*) realloc(mprimes, pc * sizeof(int));
	*primes = mprimes;
}

//Computes the Euler's totient function of n.
int phi(int n) {
	int* primes;
	int n_primes;
	allprimes(isqrt(n), &primes, &n_primes);
	int phi_n = 1;
	int p_idx = 0;
	int p;
	for (p_idx = 0; p_idx < n_primes; ++p_idx) {
		p = primes[p_idx];
		if (p > n) {
			break;
		} else {
			if (n % p == 0) {
				n /= p;
				phi_n *= p - 1;
			}

			while (n % p == 0) {
				n /= p;
				phi_n *= p;
			}
		}
	}

	if (n != 1) {
		//n must be a prime
		phi_n *= n - 1;
	}

	free(primes);
	return phi_n;
}
	
/**
 ***********************************************************
 *********************** Sequences *************************
 ***********************************************************
 */

//Returns the nth Fibonacci number in ln(n) time.
// f_0 = 1, f_1 = 1 and f_(n+2) = f_(n+1) + f_n
unsigned long fibonacci(unsigned int n) {
	if (n == 0 || n == 1) {
		return 1;
	}

	unsigned long R11 = 1L;
	unsigned long R12 = 1L;
	unsigned long R21 = 1L;
	unsigned long R22 = 0L;
	unsigned long An11 = 1L;
	unsigned long An12 = 0L;
	unsigned long An21 = 0L;
	unsigned long An22 = 1L;
	unsigned long oAn11, oAn12, oAn21, oAn22;
	unsigned long oR11, oR12, oR21, oR22;

	while (true) {
		if (n & 1) {
			oAn11 = An11;
			oAn12 = An12;
			oAn21 = An21;
			oAn22 = An22;
			An11 = oAn11 * R11 + oAn12 * R21;
			An12 = oAn11 * R12 + oAn12 * R22;
			An21 = oAn21 * R11 + oAn22 * R21;
			An22 = oAn21 * R12 + oAn22 * R22;
		}

		n >>= 1;
		if (n == 0) break;

		oR11 = R11;
		oR12 = R12;
		oR21 = R21;
		oR22 = R22;

		R11 = oR11 * oR11 + oR12 * oR21;
		R12 = oR11 * oR12 + oR12 * oR22;
		R21 = oR21 * oR11 + oR22 * oR21;
		R22 = oR21 * oR12 + oR22 * oR22;
	}
	return An11;
}

unsigned long* fiblist11(int n) {
	unsigned long* toreturn = (unsigned long *)malloc(n * sizeof(unsigned long));

	if (n == 0) {
		return toreturn;		
	}

	toreturn[0] = 1L;

	if (n == 1) {
		return toreturn;
	}

	toreturn[1] = 1L;

	if (n == 2) {
		return toreturn;
	}

	int idx;
	unsigned long fib_im1 = 1L, fib_im2 = 1L, tmp;
	for (idx = 2; idx < n; ++idx) {
		tmp = fib_im1;
		toreturn[idx] = (fib_im1 += fib_im2);
		fib_im2 = tmp;
	}

	return toreturn;
}

unsigned long* fiblist12(int n) {
	unsigned long* toreturn = (unsigned long *)malloc(n * sizeof(unsigned long));

	if (n == 0) {
		return toreturn;
	}

	toreturn[0] = 1L;

	if (n == 1) {
		return toreturn;
	}

	toreturn[1] = 2L;

	if (n == 2) {
		return toreturn;
	}

	int idx;
	unsigned long fib_im1 = 1L, fib_im2 = 2L, tmp;
	for (idx = 2; idx < n; ++idx) {
		tmp = fib_im1;
		toreturn[idx] = (fib_im1 += fib_im2);
		fib_im2 = tmp;
	}

	return toreturn;
}

unsigned int fibinv(unsigned long N) {
	const unsigned int guess = round(log(N) * LOG_GR_REC);
	if (fibonacci(guess + 1) > N) {
		return guess - 1;
	} else {
		return guess;
	}
}

/**
 ********************************************************
 ***************** Numerical schemes ********************
 ********************************************************
 */

double trapezoidal(double (*f)(double), double a, double b, int n) {
	double h = (b - a) / n;
	double I = 0.5 * ((*f)(a) + (*f)(b));
	int i;
	for (i = 1; i < n; ++i) {
		I += (*f)(a + i * h);
	}
	return I * h;
}

double rombergintegrate(double (*f)(double), double a, double b, int n, double tol) {
	double X_c0, X_l0, X_c1, error;

	X_l0 = trapezoidal(f, a, b, n);

	n <<= 1;

	X_c0 = trapezoidal(f, a, b, n);
	X_c1 = (4 * X_c0 - X_l0) / 3;

	error = fabs(X_c1 - X_l0);

	if (error < tol) {
		return X_c1;
	}

	double X_l1, X_c2;
	n <<= 1;
	X_l0 = X_c0;
	X_l1 = X_c1;

	X_c0 = trapezoidal(f, a, b, n);
	X_c1 = (4 * X_c0 - X_l0) / 3;
	X_c2 = (16 * X_c1 - X_l1) / 15;

	error = fabs(X_c2 - X_l1);

	if (error < tol) {
		return X_c2;
	}

	double X_l2, X_c3;
	n <<= 1;
	X_l0 = X_c0;
	X_l1 = X_c1;
	X_l2 = X_c2;

	X_c0 = trapezoidal(f, a, b, n);
	X_c1 = (4 * X_c0 - X_l0) / 3;
	X_c2 = (16 * X_c1 - X_l1) / 15;
	X_c3 = (64 * X_c2 - X_l2) / 63;

	error = fabs(X_c3 - X_l2);

	if (error < tol) {
		return X_c3;
	}

	double X_l3, X_c4;
	n <<= 1;
	X_l0 = X_c0;
	X_l1 = X_c1;
	X_l2 = X_c2;
	X_l3 = X_c3;

	X_c0 = trapezoidal(f, a, b, n);
	X_c1 = (4 * X_c0 - X_l0) / 3;
	X_c2 = (16 * X_c1 - X_l1) / 15;
	X_c3 = (64 * X_c2 - X_l2) / 63;
	X_c4 = (256 * X_c3 - X_l3) / 255;

	error = fabs(X_c4 - X_l3);

	if (error < tol) {
		return X_c4;
	}

	double X_l4, X_c5;
	n <<= 1;
	X_l0 = X_c0;
	X_l1 = X_c1;
	X_l2 = X_c2;
	X_l3 = X_c3;
	X_l4 = X_c4;

	X_c0 = trapezoidal(f, a, b, n);
	X_c1 = (4 * X_c0 - X_l0) / 3;
	X_c2 = (16 * X_c1 - X_l1) / 15;
	X_c3 = (64 * X_c2 - X_l2) / 63;
	X_c4 = (256 * X_c3 - X_l3) / 255;
	X_c5 = (1024 * X_c4 - X_l4) / 1023;

	error = fabs(X_c5 - X_l4);

	if (error < tol) {
		return X_c5;
	}

	double X_l5, X_c6;
	n <<= 1;
	X_l0 = X_c0;
	X_l1 = X_c1;
	X_l2 = X_c2;
	X_l3 = X_c3;
	X_l4 = X_c4;
	X_l5 = X_c5;

	X_c0 = trapezoidal(f, a, b, n);
	X_c1 = (4 * X_c0 - X_l0) / 3;
	X_c2 = (16 * X_c1 - X_l1) / 15;
	X_c3 = (64 * X_c2 - X_l2) / 63;
	X_c4 = (256 * X_c3 - X_l3) / 255;
	X_c5 = (1024 * X_c4 - X_l4) / 1023;
	X_c6 = (4096 * X_c5 - X_l5) / 4095;

	error = fabs(X_c6 - X_l5);

	if (error < tol) {
		return X_c6;
	}

	abort();
}

/**
 ******************************************************
 **************** Basic operations ********************
 ******************************************************
 */

unsigned int isqrt(unsigned int n) {
	unsigned int lower = 0;
	unsigned int upper = MIN(0xffff, n); //To avoid overflow.
	unsigned int mid;
	while (lower < upper) {
		mid = (lower + upper + 1) >> 1;
		if (mid * mid <= n) {
			lower = mid;
		} else {
			upper = mid - 1;
		}
	}
	return lower;
}

//Computes x^n where n in ln(n) time where n is a positve integer.
double dipow(double x, int n) {
	if (n == 0) {
		return 1;
	}
	double r = x;
	double y = 1;
	while (true) {
		if (n & 1) {
			y *= r;
		}
		n >>= 1;
		if (n == 0) break;
		r *= r;
	}
	return y;
}

int binsearch(int tofind, int* array, size_t len) {
	int lower = 0;
	int upper = len;
	int mid;
	while (lower < upper) {
		mid = (lower + upper) >> 1;
		if (array[mid] < tofind) {
			lower = mid + 1;
		} else {
			upper = mid;
		}
	}

	if (upper < len) {
		return array[upper] == tofind ? upper : -(upper + 1); 
	} else {
		return -(len + 1);
	}
}


// int main(int argc, char* argv[]) {
// 	/*
// 	int m,n,a,b,c;
// 	m = 899;
// 	n = 319;
// 	bezout(m, n, &a, &b, &c);
// 	printf("%d*%d + %d*%d = %d\n", a, m, b, n, c);
// 	*/
	
// 	/*
// 	int a, b;
// 	a = 2;
// 	b = 3;
// 	int a_inv = inv(a, b);
// 	printf("%d\n", a_inv);
// 	*/

// 	/*
// 	int n = 1207963139;
// 	bool ip = is_prime(n);
// 	if (ip) {
// 		printf("Success\n");
// 	} else {
// 		printf("Failure\n");
// 	}
// 	*/
	
// 	/*
// 	int* primes;
// 	int n_primes;
// 	allprimes(n, &primes, &n_primes);
// 	printf("Done\n");
// 	printf("There are %d many primes below %d \n", n_primes, n);

// 	free(primes);
// 	*/

// 	//int n = atoi(argv[1]);
// 	//unsigned long f_n = fibonacci(n);
// 	//printf("fibonacci(%d) = %lu\n", n, f_n);

	
// 	int idx;
// 	int tofind;

// 	const int len = 14;
// 	int arr1[] = {0, 0, 0, 1, 2, 3, 4, 4, 6, 7, 8, 9, 9, 9};
// 	tofind = 0;
// 	idx = binsearch(tofind, arr1, len);
// 	printf("Index of %d is %d\n", tofind, idx);

// 	tofind = 9;
// 	idx = binsearch(tofind, arr1, len);
// 	printf("Index of %d is %d\n", tofind, idx);

// 	tofind = 4;
// 	idx = binsearch(tofind, arr1, len);
// 	printf("Index of %d is %d\n", tofind, idx);

// 	tofind = 6;
// 	idx = binsearch(tofind, arr1, len);
// 	printf("Index of %d is %d\n", tofind, idx);

// 	tofind = 10;
// 	idx = binsearch(tofind, arr1, len);
// 	printf("Index of %d is %d\n", tofind, idx);

// 	tofind = -1;
// 	idx = binsearch(tofind, arr1, len);
// 	printf("Index of %d is %d\n", tofind, idx);

// 	int n = atoi(argv[1]);
// 	unsigned int k = fibinv(n);
// 	unsigned long fibk = fibonacci(k + 1);
// 	printf("%d\n", n);
// 	printf("%d\n", k);
// 	printf("%lu\n", fibk);

// }
