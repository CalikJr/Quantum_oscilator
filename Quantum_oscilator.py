import numpy 
import matplotlib.pyplot as plt
def osc(x):
    return x**2

x_min = float(input())
x_max = float(input())
no = int(input())

step = (x_max-x_min)/no

Grid = numpy.linspace(x_min,x_max,no)

H = numpy.zeros([no,no])

for i in range(no):
        H[i, i] = (2.0 /(step * step) + (2.0 * osc(Grid[i]) * 1 / (1 * 1))  )
        if i > 0:
            H[i, i - 1] = -1.0/ (step * step)
        if i < (no - 1):
            H[i, i + 1] = -1.0  / (step * step)


E, vects = numpy.linalg.eig(H)
print(E)
print(vects)
plt.plot(Grid,vects[:,0])
plt.show()