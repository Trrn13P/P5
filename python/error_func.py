import numpy as np
import matplotlib.pyplot as plt

directory = "../textfiles/"

filenames = []
filenames.append(directory+"forward_euler.txt")
filenames.append(directory+"backward_euler.txt")
filenames.append(directory+"crank_nicolson.txt")


def u_analytic(t_,m_max,n,x):
    u_analytic_ = np.zeros(n+1)
    for m in range(1,m_max):
        A_m = 2/(m*np.pi)*(1-np.cos(m*np.pi))
        u_analytic_ += A_m*np.sin(m*np.pi*x)*np.exp(-m**2*np.pi**2*t_)
    return u_analytic_

u_all = []
schemes = []

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
    t_val = dt*(saved_tsteps+1)
    t = np.linspace(0,t_val,len(u))
    scheme = filename[13:].split(".")[0].split("_")[0]

    u_all.append(u)
    schemes.append(scheme)


epsilon_forward = []
epsilon_backward = []
epsilon_cranck = []
m_max = 100

epsilons = [epsilon_forward,epsilon_backward,epsilon_cranck]

for u_scheme,k in zip(u_all,[0,1,2]):
    for j in range(0,len(t)):
        local_epsilon = 0
        u_analytic_ = u_analytic(t[j],m_max,n,x)
        for i in range(0,n+1):
            local_epsilon += abs(u_analytic_[i]-u_scheme[j][i])/(n+1)
        epsilons[k].append(local_epsilon)

plt.plot(t,epsilons[0],label="forward")
plt.plot(t,epsilons[1],label="backward")
plt.plot(t,epsilons[2],label="crank")
plt.xlabel("t")
plt.ylabel(r"$\varepsilon$(t)")
plt.legend()
#plt.show()

dir = "../figures/d/"
filename_ = "n=" + str(n) + ",dt=" + str(dt) +",dx=" +str(dx)+ ", epsilon, M=" +str(m_max) + ",T="+str(t[-1])
plt.savefig(dir+filename_+".png")
