import numpy as np
import matplotlib.pyplot as plt

directory = "../textfiles/"



#creating array for all filenames made from c++ program, only for cranc-nicolson
filenames = []
filenames.append(directory+"crank_nicolson.txt")


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
    scheme = filename[13:].split(".")[0].split("_")[0]

#plotting crank_nicolson for the 50-th timestep, also doing this for the
#analytic solutions
for index in [50]:
    t_val = index*dt*(saved_tsteps+1)
    plt.plot(x,u[index],label=scheme)

    #creating analytic solution for different m_max, which is the
    #max of the sum for the analytic solution
    for m_max in [2,11,10001]:
        u_analytic = np.zeros(n+1)
        for m in range(1,m_max):
            A_m = 2/(m*np.pi)*(1-np.cos(m*np.pi))
            u_analytic += A_m*np.sin(m*np.pi*x)*np.exp(-m**2*np.pi**2*t_val)
        plt.plot(x,u_analytic,label=r"analytic $M$=%g"%(m_max-1))
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
#saving file
dir = "../figures/d/"
filename_ = "n=" + str(n) + ",dt=" + str(dt) +",dx=" +str(dx)+ ",t=" + str(t_val) + ",crank"
plt.savefig(dir+filename_+".png")
plt.clf()
