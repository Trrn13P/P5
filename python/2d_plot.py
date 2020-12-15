import numpy as np
import matplotlib.pyplot as plt

directory = "../textfiles/"

filenames = []
filenames.append(directory+"forward_euler.txt")
filenames.append(directory+"backward_euler.txt")
filenames.append(directory+"crank_nicolson.txt")


"""

for filename in filenames:
    infile = open(filename,"r")
    first_line = infile.readline().split()
    dt = eval(first_line[0].split("=")[1])
    saved_tsteps = eval(first_line[1].split("=")[1])
    alpha = eval(first_line[2].split("=")[1])
    dx = eval(first_line[3].split("=")[1])
    infile.readline()

    u = []
    for line in infile:
        u.append([eval(element) for element in line.split()])
    infile.close()


    n = len(u[0])-1
    x = np.linspace(0,1,n+1)

    scheme = filename[13:].split(".")[0].split("_")[0]
    for index in [5,10]:
        plt.plot(x,u[index],label=scheme+", t="+str(index*dt*saved_tsteps) + "s")

plt.legend(fontsize=10)
savefig_filename = "n=" +str(n)+",dt="+str(dt)+",dx="+str(dx)
dir = "../figures/c/"
plt.savefig(dir+savefig_filename+".png")
plt.clf()
#plt.show()


"""


#t = np.linspace(0,dt*(n+1),n)
n= 1000
x = np.linspace(0,1,n)
dt = 0.001

for t in [0,100*dt,200*dt]:
    u = np.zeros(n)
    for m in range(1,10000):
        A_m = 2/(m*np.pi)*(1-np.cos(m*np.pi))
        u += A_m*np.sin(m*np.pi*x)*np.exp(-m**2*np.pi**2*t)
    plt.plot(x,u,label="t=%g"%(t))
plt.legend()
plt.show()
