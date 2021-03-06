import math


############################################################
# Various utility functions for mathematical problems. 
# Author: Hjalti Thor Isleifsson
############################################################

############################################################
#################### Number theory #########################
############################################################

def gcd(a, b):
	while b != 0:
		q = a // b
		r = a - q * b
		a = b
		b = r

	return a

def lcm(a, b):
	return a * b // gcd(a, b)

#Returns numbers a,b such that a * m + b * n = gcd(m,n)
def bezout(m, n):
	q = m // n
	old_r = n
	r = m - q * n

	old_a = 1
	a = 0
	old_b = 0
	b = 1

	while r != 0:
		old_old_a = old_a
		old_a = a
		a = old_old_a - q * a

		old_old_b = old_b
		old_b = b
		b = old_old_b - q * b

		q = old_r // r
		
		old_old_r = old_r
		old_r = r
		r = old_old_r - q * r

	return (a, b, old_r)

def inv(x, p):
	(a, b, gcd) = bezout(x, p)
	if gcd != 1:
		raise Exception('%d does not have an inverse in Z_%d' % (x,p))
	elif a < 0:
		return a + p
	else:
	 	return a

def isPrime(n):
	if n < 4: 
		return n > 1
	elif n % 2 == 0 or n % 3 == 0:
		return False
	else:
		d = 5
		while d * d <= n:
			if n % d == 0 or n % (d + 2) == 0:
				return False

			d += 6

		return True


#Returns a list of all primes less then or equal to n.
def allPrimes(n):
	if n < 2:
		return []
	if n < 3:
		return [2]
	if n < 5: 
		return [2, 3]

	#We use segmented sieve.
	primes = [0] * int(1.2 * n / math.log(n))
	p_idx = 0
	d = sqrt(n) #Length of the segment.
	flags = [True] * d
	flags[0] = False
	flags[1] = False

	for i in range(0, d):
		if flags[i]:
			primes.insert(p_idx, i)
			p_idx += 1
			for j in range(2 * i, d, i):
				flags[j] = False

	n_iv = (n + d - 1) // d #Number of intervals
	for iv_idx in range(1, n_iv):
		offset = iv_idx * d
		upto = min(n + 1, offset + d)
		sqrt_upto = sqrt(upto)

		for i in range(d):
			flags[i] = True

		for p in primes:
			if p > sqrt_upto or p == 0:
				break
			else:
				for j in range(((offset + p - 1) // p) * p - offset, upto - offset, p):
					flags[j] = False

		for i in range(upto - offset):
			if flags[i]:
				primes.insert(p_idx, i + offset)
				p_idx += 1

	return primes[0:p_idx]

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

#######################################################
################# Basic operations ####################
#######################################################

def sqrt(n):
	if n == 0:
		return 0
	elif n < 4:
		return 1
	elif n < 9:
		return 2
	elif n < 16:
		return 3
	elif n < 25:
		return 4
	else:
		lower = 5
		upper = n // 5
		while lower < upper: 
			x = (lower + upper + 1) >> 1
			if x * x <= n:
				lower = x
			else:
				upper = x - 1

		return lower

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
