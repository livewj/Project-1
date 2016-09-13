#plots the results produced in "example.cpp"
#where n=10, n=100 and n=1000

import numpy as np
import matplotlib.pyplot as plt

#function that reads the data produced in "example.cpp"

def read(filename):
	infile = open(filename, 'r')
	x = []
	u = []
	v = []
	relevant_lines = infile.readlines()[1:] #skips the first line
	for line in relevant_lines:
		data = line.split()
		x.append(float(data[0]))
		u.append(float(data[1]))
		v.append(float(data[2]))
	infile.close()
	x = np.array(x)
	u = np.array(u)
	v = np.array(v)
	return x, u, v

#n=10
x1, u1, v1 = read('data_n10')
#n=100
x2, u2, v2 = read('data_n100')
#n=1000
x3, u3, v3 = read('data_n1000')

#plotting
plt.plot(x1,u1, 'r--', x1, v1, 'b')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend(['Analytical', 'Numerical'])
plt.title('Numerical and analytical solution for n = 10')
plt.savefig('b_n10.png')
plt.show()

plt.plot(x2,u2, 'r--', x2, v2, 'b')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend(['Analytical', 'Numerical'])
plt.title('Numerical and analytical solution for n = 100')
plt.savefig('b_n100.png')
plt.show()

plt.plot(x3,u3, 'r--', x3, v3, 'b')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend(['Analytical', 'Numerical'])
plt.title('Numerical and analytical solution for n = 1000')
plt.savefig('b_n1000.png')
plt.show()
