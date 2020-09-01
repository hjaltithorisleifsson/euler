import math
import math_utils as mu

def main():
	a = 3 / 5
	b = -12 / 5
	c = math.pi - 2
	d = b**2 - 4 * a * c
	v = (-b - math.sqrt(d)) / (2 * a)

	bv_prime_norm = lambda t: B_prime_norm(t, v)
	tol = 0.5 * (10**(-10))
	L = mu.integrate(bv_prime_norm, 0, 1, tol)

	r = 100 * (L - 0.5 * math.pi) / (0.5 * math.pi)
	print(r)
	

def B_prime_norm(t,v):
	return math.sqrt(9*((1-t)**4) + 9*((1-t)**2)*((1-3*t)**2)*(1+v*v) + 9*((2-3*t)**2)*(t**2)*(1+v*v) + 9*(t**4) - 18*((1-t)**3)*(1-3*t)-18*((1-t)**2)*(2-3*t)*t*v + 36*(1-t)*(1-3*t)*(2-3*t)*t*v+18*(1-t)*(1-3*t)*(t**2)*v + 18*(2-3*t)*(t**3))

main()