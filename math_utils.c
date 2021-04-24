#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (y)) ? (x) : (y))

/*
############################################################
# Various utility functions for mathematical problems. 
# Author: Hjalti Thor Isleifsson
############################################################

############################################################
#################### Number theory #########################
############################################################
*/

int gcd(int a, int b);
	

int lcm(int a, int b);


void bezout(int m, int n, int* a, int* b, int* c);


int inv(int x, int p);


bool is_prime(int n);

//All primes <= n
void all_primes(int n, int** primes, int* n_primes);

//Computes the Euler's totient function of n.
int phi(int n);

/*
############################################################
######################## Sequences #########################
############################################################
*/


//Returns the nth Fibonacci number in ln(n) time.
//f_0 = 1, f_1 = 1 and f_(n+2) = f_(n+1) + f_n
long fibonacci(int n);

/*
#f: An increasing unbounded sequence.
#Finds the greatest i s.t. f_i <= n using binary search.
def bin_inv_seq(f, n):
	lower = 0
	upper = 1
	while f(upper) <= n:
		lower = upper
		upper <<= 1

	while lower < upper: 
		mid = (lower + upper + 1) >> 1
		if f(mid) <= n:
			lower = mid
		else:
			upper = mid - 1

	return lower*/
/*
#########################################################
################## Numerical schemes ####################
#########################################################
*/
/*
def trapezoidal(f, a, b, n):
	h = (b - a) / n
	I = 0.5 * (f(a) + f(b))
	for i in range(1, n):
		I += f(a + i * h)

	return I * h
	

//double integrate(f, a, b, tol); 
*/

/*
#######################################################
################# Basic operations ####################
#######################################################
*/

unsigned int isqrt(unsigned int n);

//Computes x^n where n in ln(n) time where n is a positve integer.
int dipow(double x, int n);

//#####################################################
// Here come the implementations
//#####################################################

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
		return(1);
	} else if (a < 0) {
		return a + p;
	} else {
	 	return a;
	}
}

bool is_prime(int n) {
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
void all_primes(int n, int** primes, int* n_primes) {
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
	mprimes = (int *)realloc(mprimes, pc * sizeof(int));
	*primes = mprimes;
}

/*
#Computes the Euler's totient function of n.
def phi(n):
	sqrt_n = sqrt(n)
	primes = allPrimes(sqrt_n)
	phi_n = 1

	for p in primes:
		if p > n: 
			break
		else:
			if n % p == 0:
				n //= p
				phi_n *= p - 1

			while n % p == 0:
				n //= p
				phi_n *= p

	if n != 1:
		#n must be a prime
		phi_n *= n - 1

	return phi_n

############################################################
######################## Sequences #########################
############################################################


#Returns the nth Fibonacci number in ln(n) time.
# f_0 = 1, f_1 = 1 and f_(n+2) = f_(n+1) + f_n
def fibonacci(n):
	if n == 0 or n == 1:
		return 1

	R_11 = 1
	R_12 = 1
	R_21 = 1
	R_22 = 0
	A_n_11 = 1
	A_n_12 = 0
	A_n_21 = 0
	A_n_22 = 1

	while True:
		if n & 1 == 1:
			old_A_n_11 = A_n_11
			old_A_n_12 = A_n_12
			old_A_n_21 = A_n_21
			old_A_n_22 = A_n_22
			A_n_11 = old_A_n_11 * R_11 + old_A_n_12 * R_21
			A_n_12 = old_A_n_11 * R_12 + old_A_n_12 * R_22
			A_n_21 = old_A_n_21 * R_11 + old_A_n_22 * R_21
			A_n_22 = old_A_n_21 * R_12 + old_A_n_22 * R_22

		n >>= 1
		if n == 0:
			break

		old_R_11 = R_11
		old_R_12 = R_12
		old_R_21 = R_21
		old_R_22 = R_22

		R_11 = old_R_11 * old_R_11 + old_R_12 * old_R_21
		R_12 = old_R_11 * old_R_12 + old_R_12 * old_R_22
		R_21 = old_R_21 * old_R_11 + old_R_22 * old_R_21
		R_22 = old_R_21 * old_R_12 + old_R_22 * old_R_22

	return A_n_11

#f: An increasing unbounded sequence.
#Finds the greatest i s.t. f_i <= n using binary search.
def bin_inv_seq(f, n):
	lower = 0
	upper = 1
	while f(upper) <= n:
		lower = upper
		upper <<= 1

	while lower < upper: 
		mid = (lower + upper + 1) >> 1
		if f(mid) <= n:
			lower = mid
		else:
			upper = mid - 1

	return lower

#########################################################
################## Numerical schemes ####################
#########################################################

def trapezoidal(f, a, b, n):
	h = (b - a) / n
	I = 0.5 * (f(a) + f(b))
	for i in range(1, n):
		I += f(a + i * h)

	return I * h

def integrate(f, a, b, tol):
	error = math.inf
	max_it = 100
	it = 0

	n = 2
	X_i = [trapezoidal(f,a,b,n)]
	while error > tol and it < max_it:
		i = len(X_i)
		n *= 2
		X_ip = [0 for i in range(i+1)]
		X_ip[0] = trapezoidal(f,a,b,n)
		for j in range(1, i + 1):
			fj = 4**j
			X_ip[j] = (fj * X_ip[j-1] - X_i[j-1]) / (fj - 1)

		error = math.fabs(X_i[-1] - X_ip[-1])
		it += 1
		X_i = X_ip

	return X_i[-1]
*/

/*
#######################################################
################# Basic operations ####################
#######################################################
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

/*
#Computes x^n where n in ln(n) time where n is a positve integer.
def pow(x, n):
	if n == 0:
		return 1

	r = x
	y = 1
	while True:
		if n & 1 == 1:
			y *= r
		
		n >>= 1

		if n == 0:
			break

		r *= r

	return y
*/

int main(int argc, char* argv[]) {
	/*
	int m,n,a,b,c;
	m = 899;
	n = 319;
	bezout(m, n, &a, &b, &c);
	printf("%d*%d + %d*%d = %d\n", a, m, b, n, c);
	*/
	
	/*
	int a, b;
	a = 2;
	b = 3;
	int a_inv = inv(a, b);
	printf("%d\n", a_inv);
	*/

	/*
	int n = 1207963139;
	bool ip = is_prime(n);
	if (ip) {
		printf("Success\n");
	} else {
		printf("Failure\n");
	}
	*/
	/*
	int n = atoi(argv[1]);
	int* primes;
	int n_primes;
	all_primes(n, &primes, &n_primes);

	int i;
	for (i = 0; i < n_primes; ++i) {
		printf("%d\n", primes[i]);
	}

	free(primes);*/
}