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


#Returns a list of all primes less then or equal to n.
def allPrimes(n):
	primes = [0] * int(1.2 * n / math.log(n))
	p_idx = 0
	d = sqrt(n)
	flags = [True] * d
	flags[0] = False
	flags[1] = False

	for i in range(2, d):
		if flags[i]:
			primes.insert(p_idx, i)
			p_idx += 1
			for j in range(2 * i, d, i):
				flags[j] = False
		
		else:
			flags[i] = True

	n_iv = n // d + 1
	for iv_idx in range(1, n_iv):
		offset = iv_idx * d
		upto = min(n + 1, offset + d)

		sqrt_upto = sqrt(upto)
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
			else:
				flags[i] = True

	return primes[0:p_idx]

#Returns the nth Fibonacci number
# f_0 = 1, f_1 = 1 and f_(n+2) = f_(n+1) + f_n
def fibonacci(n):
	if n == 0 or n == 1:
		return 1
	else: 
		A_11 = 1
		A_12 = 1
		A_21 = 1
		A_22 = 0

		A_n_11 = 1
		A_n_12 = 1
		A_n_21 = 1
		A_n_22 = 0

		m = n - 2
		while m != 0:
			A_nm_11 = A_n_11
			A_nm_12 = A_n_12
			A_nm_21 = A_n_21
			A_nm_22 = A_n_22

			if m & 1 == 1:
				A_n_11 = A_nm_11 * A_11 + A_nm_12 * A_21
				A_n_12 = A_nm_11 * A_12 + A_nm_12 * A_22
				A_n_21 = A_nm_21 * A_11 + A_nm_22 * A_21
				A_n_22 = A_nm_21 * A_12 + A_nm_22 * A_22
				m -= 1
			else:
				A_n_11 = A_nm_11 * A_nm_11 + A_nm_12 * A_nm_21
				A_n_12 = A_nm_11 * A_nm_12 + A_nm_12 * A_nm_22
				A_n_21 = A_nm_21 * A_nm_11 + A_nm_22 * A_nm_21
				A_n_22 = A_nm_21 * A_nm_12 + A_nm_22 * A_nm_22
				m >>= 1

		return A_n_11 + A_n_12




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
