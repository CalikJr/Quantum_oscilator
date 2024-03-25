import numpy 
import matplotlib.pyplot as plt
def osc(x):
    return x**2

x_min = -15#float(input())
x_max = 15#float(input())
no = 5500#int(input())
rm = 4

step = (x_max-x_min)/no

Grid = numpy.linspace(x_min,x_max,no)

H = numpy.zeros([no,no])

for i in range(no):
        H[i, i] = ((2.0 /(step * step) +  osc(Grid[i]))) 
        if i > 0:
            H[i, i - 1] = -1.0/ (step * step)
        if i < (no - 1):
            H[i, i + 1] = -1.0  / (step * step)


E, vects = numpy.linalg.eigh(H)
#E /= (step*step)
E /= 2
print(E[rm])
#print(vects)
plt.plot(Grid,vects[:,rm])
plt.show()