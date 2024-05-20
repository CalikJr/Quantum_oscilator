import numpy 
import matplotlib.pyplot as plt
from matplotlib import animation
import math
import sympy

def osc(x, epsilon):
    return x**2+epsilon*x**4

x_min = -15#float(input())
x_max = 15#float(input())
no = 500#int(input())
rm = 4

step = (x_max-x_min)/no

Grid = numpy.linspace(x_min,x_max,no)

def energy(Epsilon):
    H = numpy.zeros([no,no])
    for i in range(no):
            H[i, i] = ((2.0 /(step * step) +  osc(Grid[i], Epsilon))) 
            if i > 0:
                H[i, i - 1] = -1.0/ (step * step)
            if i < (no - 1):
                H[i, i + 1] = -1.0  / (step * step)
    E, vects = numpy.linalg.eigh(H)
    #E /= (step*step)
    E /= 2
    #print(E[rm])
    #print(vects)
    #plt.plot(Grid,vects[:,rm])
    #plt.show()
    vects=vects
    return E, vects

E, vects = energy(0)

C0=1

def psi0(x):
    return C0*math.exp(-(x**2)/2)



def H(x, n):
    match n:
        case 0:
            return 1
        case 1:
            return 2*x
        case 2:
            return 4*(x**2)-2
        case 3:
            return 8*(x**3)-12*x
        case 4:
            return 16*(x**4)-48*(x**2)+12
        case 5:
            return 32*(x**5)-160*(x**3)+120*x

Cn = 1
def psin(x, n):
    return Cn*H(x, n)*psi0(x)

for i in range(6):
    Hvect=numpy.zeros(Grid.size)
    gn = 0
    sum = 0
    for g in Grid:
        sum += psin(g, i)*psin(g, i)
    Cn = 1.0/math.sqrt(sum)
    print(Cn)
    for g in Grid:
        Hvect[gn]=psin(g, i)
        gn+=1
    Cn=1
    plt.plot(Grid, vects[:, i], label="Obliczenia")
    plt.plot(Grid, Hvect, 'r.', ms=3, label="Wielomian")
    plt.legend(loc="upper left")
    plt.title("Wektor własny stanu podstawowego E"+str(i))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

epsilon_range = numpy.linspace(0, 2, 100)
E0 = numpy.zeros([100])
E1 = numpy.zeros([100])
E2 = numpy.zeros([100])
E3 = numpy.zeros([100])
n=0
for i in epsilon_range:
    E, vects = energy(i)
    #print(E[0])
    E0[n] = E[0]
    E1[n] = E[1]
    E2[n] = E[2]
    E3[n] = E[3]
    n+=1
plt.plot(epsilon_range, E0, label="E0")
plt.plot(epsilon_range, E1, label="E1")
plt.plot(epsilon_range, E2, label="E2")
plt.plot(epsilon_range, E3, label="E3")
plt.legend(loc="upper left")
plt.title("Wykres wartości włąsnych od epsilona")
plt.xlabel("epsilon")
plt.ylabel("wartość energii")
plt.show()
print(E1)
#E, vects = energy(0)
a=0.5
sigma=1
P=numpy.array([])
for x in Grid:
    P=numpy.append(P, (1/(sigma*math.sqrt(2*math.pi)))*math.exp(-0.5*((x-a)/sigma)**2))
print(vects.shape)
print(P.size)
#plt.plot(Grid, P)
#plt.show()
c_const=numpy.array([])
for i in range(no):
    c_const=numpy.append(c_const, step*numpy.dot(P, vects[:,i]))
#print(c_const.size)
#print(c_const)

def packet_state(t):
    Packet=numpy.array([])
    for i in range(no):
        sum=0
        for j in range(no):
            sum+=c_const[j]*(math.cos(E[j]*t)+1j*math.sin(E[j]*t))*vects[i, j]
        Packet=numpy.append(Packet, abs(sum))
    #print(Packet.size)
    #print(Packet.shape)
    return Packet
#plt.plot(Grid, Packet[:])
#plt.show()

t_final = 10.0
t_initial = 0.0
t = t_initial
dt = 1

fig, ax = plt.subplots()
ax.set_xlabel('x')
ax.set_ylabel('P')
plotLine, = ax.plot(Grid, numpy.zeros(len(Grid))*numpy.NaN, 'r-')
plotTitle = ax.set_title("t=0")
ax.set_ylim(-0.01,0.01)
ax.set_xlim(-4,4)

def animate(t):
    pp = packet_state(t)
    plotLine.set_ydata(pp)
    plotTitle.set_text(f"t = {t:.1f}")
    #ax.relim() # use if autoscale desired
    #ax.autoscale()
    return [plotLine,plotTitle]



ani = animation.FuncAnimation(fig, func=animate, frames=numpy.arange(t_initial, t_final+dt, dt), blit=True)
plt.plot(Grid, osc(Grid, 0))
plt.show()