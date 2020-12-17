import numpy as np
import matplotlib.pyplot as plt

directory = "../textfiles/"

#creating array for all filenames made from c++ program
filenames = []
filenames.append(directory+"forward_euler.txt")
filenames.append(directory+"backward_euler.txt")
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

    #plotting for different t-values
    scheme = filename[13:].split(".")[0].split("_")[0]
    timestep_numbers = [5,10]
    for index in timestep_numbers:
        plt.plot(x,u[index],label=scheme+", t="+str(index*dt*(saved_tsteps+1)) + "s")

#saving file
plt.legend(fontsize=10)
savefig_filename = "n=" +str(n)+",dt="+str(dt)+",dx="+str(dx)
dir = "../figures/c/"
plt.savefig(dir+savefig_filename+".png")
plt.clf()
