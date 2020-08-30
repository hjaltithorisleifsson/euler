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
