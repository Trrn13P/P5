import numpy as np
import matplotlib.pyplot as plt

directory = "../textfiles/"

#creating array for all filenames made from c++ program
filenames = []
filenames.append(directory+"forward_euler.txt")
filenames.append(directory+"backward_euler.txt")
filenames.append(directory+"crank_nicolson.txt")


#analytic solution for timestep t_, M=m_max, xsteps n and x-array
def u_analytic(t_,m_max,n,x):
    u_analytic_ = np.zeros(n+1)
    for m in range(1,m_max):
        A_m = 2/(m*np.pi)*(1-np.cos(m*np.pi))
        u_analytic_ += A_m*np.sin(m*np.pi*x)*np.exp(-m**2*np.pi**2*t_)
    return u_analytic_

#u_all will have 3 elements, first is the u-values for all timesteps of forward,
#then all for backward, and then all for crank-nicolson
u_all = []
schemes = []

for filename in filenames:
    #opening infile and saving elements on first line
    infile = open(filename,"r")
    first_line = infile.readline().split()
    dt = eval(first_line[0].split("=")[1])
    saved_tsteps = eval(first_line[1].split("=")[1])
    alpha = eval(first_line[2].split("=")[1])
    dx = eval(first_line[3].split("=")[1])
    infile.readline()

    #appending u-values for each of the timesteps, every line is a new timestep
    u = []
    for line in infile:
        u.append([eval(element) for element in line.split()])
    infile.close()

    #setting up array for x
    n = len(u[0])-1
    x = np.linspace(0,1,n+1)
    t_val = dt*(saved_tsteps+1)
    t = np.linspace(0,t_val,len(u))
    scheme = filename[13:].split(".")[0].split("_")[0]

    #appinding all timesteps u-values
    u_all.append(u)
    schemes.append(scheme)

#creating epsilon arrays for all schemes
epsilon_forward = []
epsilon_backward = []
epsilon_cranck = []
m_max = 100

#making list to hold epsilon values
epsilons = [epsilon_forward,epsilon_backward,epsilon_cranck]

#going through the u-values in u_all,  u(x,t_j) and creating a local epsilon at that timestep,
#and then summing over all u(x_i,t_j) (i's) at the timestep to find the total error.
#dividing by len(u(x,t))=(n+1) to find average error
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

#saving to png
dir = "../figures/d/"
filename_ = "n=" + str(n) + ",dt=" + str(dt) +",dx=" +str(dx)+ ", epsilon, M=" +str(m_max) + ",T="+str(t[-1])
#plt.savefig(dir+filename_+".png")
plt.show()
