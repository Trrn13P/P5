import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def plot_3d(x,t,u,save_filename):
    X,T = np.meshgrid(x,t)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, T, u, rstride=1, cstride=1, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False,shade=False)
    surf.set_facecolor((0,0,0,0))
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("u(x)")
    ax.view_init(elev=20., azim=40)

    dir = "../figures/"
    save_filename = dir + save_filename + ".png"
    plt.savefig(save_filename)
    plt.clf()






directory = "../textfiles/"

filenames = []
filenames.append(directory+"forward_euler.txt")
filenames.append(directory+"backward_euler.txt")
filenames.append(directory+"crank_nicolson.txt")




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


    x = np.asarray(x)
    u = np.asarray(u)
    t = np.linspace(0,dt*(saved_tsteps+1)*len(u),len(u))

    scheme = filename[13:].split(".")[0].split("_")[0]
    save_filename = "n="+str(n)+",dt="+str(dt)+",dx="+str(dx) + ","+scheme
    plot_3d(x,t,u,save_filename)
