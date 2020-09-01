import math

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
		X_i = X_ip

	return X_i[-1]
