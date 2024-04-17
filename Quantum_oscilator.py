import numpy 
import matplotlib.pyplot as plt
from matplotlib import animation
import math

def osc(x):
    return x**2

x_min = -15#float(input())
x_max = 15#float(input())
no = 1000#int(input())
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
#plt.plot(Grid,vects[:,rm])
#plt.show()
vects=vects/2

a=3
sigma=3
P=numpy.array([])
for x in Grid:
    P=numpy.append(P, (1/(sigma*math.sqrt(2*math.pi)))*math.exp(-0.5*((x-a)/sigma)**2))
print(vects.shape)
print(P.size)
plt.plot(Grid, P)
plt.show()
c_const=numpy.array([])
for i in range(no):
    c_const=numpy.append(c_const, numpy.dot(P, vects[:,i]))
#print(c_const.size)
#print(c_const)

def packet_state(t):
    Packet=numpy.array([])
    for i in range(no):
        sum=0
        for j in range(no):
            sum+=c_const[j]*(math.cos(E[j]*t)+math.sin(E[j]*t))*vects[i, j]
        Packet=numpy.append(Packet, sum)
    #print(Packet.size)
    #print(Packet.shape)
    return Packet
#plt.plot(Grid, Packet[:])
#plt.show()

t_final = 10.0
t_initial = 0.0
t = t_initial
dt = 0.1

fig, ax = plt.subplots()
ax.set_xlabel('x')
ax.set_ylabel('P')
plotLine, = ax.plot(Grid, numpy.zeros(len(Grid))*numpy.NaN, 'r-')
plotTitle = ax.set_title("t=0")
ax.set_ylim(-0.3,0.3)
ax.set_xlim(x_min,x_max)

def animate(t):
    pp = packet_state(t)
    plotLine.set_ydata(pp)
    plotTitle.set_text(f"t = {t:.1f}")
    #ax.relim() # use if autoscale desired
    #ax.autoscale()
    return [plotLine,plotTitle]



ani = animation.FuncAnimation(fig, func=animate, frames=numpy.arange(t_initial, t_final+dt, dt), blit=True)
plt.show()