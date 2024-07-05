# generates the table for L+1 bits of reduction. The table has size
# sizeof(dd) * 2^L
def generate_table(L):
	R = RealField(256)
	R53 = RealField(53)
	print("static const")
	print("double t[" + str(2**L) + "][2] = {")
	for i in range(2**L):
		x = R(log2(1 + i/2**L))
		a = R53(x)
		b = R53(x - a.exact_rational())
		print("\t{"+get_hex(a)+","+get_hex(b)+"},")
	print("};")
