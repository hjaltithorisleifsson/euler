import math
import math_utils as mu

def main():
	a = 1504170715041707
	b = 4503599627370517
	mod = 0
	c_sum = 0
	c = b
	for n in range(1, 50000000):
		mod = (mod + a) % b
		if mod < c:
			c = mod
			c_sum += mod

	upto = c
	print(upto)
	a_rec = mu.inv(a, b)
	ln = a_rec
	i_rec = a_rec
	c_sum += 1 
	for i in range(2, upto):
		i_rec = i_rec + a_rec
		if i_rec > b:
			i_rec -= b

		if i_rec < ln:
			c_sum += i
			ln = i_rec

	print(c_sum)

main()