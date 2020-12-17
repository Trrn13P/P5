import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

#creates 3d plot with x,t,u on axis
def plot_3d(x,t,u,save_filename):
    X,T = np.meshgrid(x,t)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, T, u, rstride=1, cstride=1, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False,shade=False)
    surf.set_facecolor((0,0,0,0))
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.set_xlim3d(0, 1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(0,1)
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("u(x)")
    ax.view_init(elev=20., azim=40)

    dir = "../figures/e/"
    save_filename = dir + save_filename + ".png"
    plt.savefig(save_filename)
    plt.clf()



# G and F as described in theory
def F(x,y,n,m):
    return np.sin(m*np.pi*x)*np.sin(n*np.pi*y)

def G(t,lmd):
    return np.exp(-lmd**2 * t)


#returning analytic expression at one timestep
def u_analytic(x,y,m_max, t_, n_len):
    u_analytic = np.zeros((n_len+1,n_len+1))
    n_max = m_max
    for i in range(0,n_len+1):
        for j in range(0,n_len+1):
            for m in range(1,m_max):
                for n in range(1,n_max):
                    lmd = np.pi * np.sqrt(m**2+n**2)
                    A_nm = 4/(n*m*np.pi**2)*(1-np.cos(m*np.pi))*(1-np.cos(n*np.pi))
                    u_analytic[i,j] += A_nm * F(x[i],y[j],n,m) * G(t_, lmd)

    return u_analytic



#opening infile
infile = open("../textfiles/2d_euler.txt")
first_line = infile.readline()
infile.readline()

#saving values from first line
dt = eval(first_line.split(" ")[0].split("=")[1])
saved_tsteps = eval(first_line.split(" ")[1].split("=")[1])
alpha = eval(first_line.split(" ")[2].split("=")[1])
dx = eval(first_line.split(" ")[3].split("=")[1])
n = eval(first_line.split(" ")[4].split("=")[1])

#creating list for all timesteps, and array for one timestep
u_all = []
u = np.zeros((n+1,n+1))

#saving all values for u, for each timestep from infile
i = 0
for line in infile:
    line = [eval(element) for element in line.split()]
    if line == []:
        u_all.append(u)
        u = np.zeros((n+1,n+1))
        i = 0
    else:
        u[i,:]=line
        i+=1
infile.close()

#position arrays
x = np.linspace(0,1,n+1)
y = np.linspace(0,1,n+1)


#max value the analytic sum will go to
m_max = 100
analytical_values = []

#time array
t = np.linspace(0,len(u_all)*dt*(saved_tsteps+1),len(u_all))

#making analytic values
for i in range(0,len(u_all)):
    analytical_values.append(u_analytic(x,y,m_max, t[i],n))

#creating the error epsilon(t)
epsilon = []
for k in range(0,len(u_all)):
    local_epsilon = 0
    for i in range(n+1):
        for j in range(n+1):
            local_epsilon += abs(analytical_values[k][i][j]-u_all[k][j][i])/((n+1)**2)
    epsilon.append(local_epsilon)
plt.plot(t,epsilon)
plt.xlabel("t")
plt.ylabel(r"$\varepsilon$(t)")

dir = "../figures/e/"
filename = "epsilon(t),n=" + str(n) +",dt=" +str(dt) + ",dx=" +str(dx) + ",M=" + str(m_max)
plt.savefig(dir+filename+".png")
plt.clf()


#making 3d plots of the analytical and numerical solutions at a couple of timesteps
for i in [0,int(len(t)/8),int(len(t)/4)]:
    filename = "3d_analytic,n=" + str(n) +",dt=" +str(dt) + ",dx=" +str(dx) + ",M=" + str(m_max) + ",T=" +str(t[i])
    plot_3d(x,y,analytical_values[i],filename)

for i in [0,int(len(t)/8),int(len(t)/4)]:
    filename = "3d_numeric,n=" + str(n) +",dt=" +str(dt) + ",dx=" +str(dx) + ",M=" + str(m_max) + ",T=" +str(t[i])
    plot_3d(x,y,u_all[i],filename)

#meshgrid for contour-plot
X,Y = np.meshgrid(x,y)

#creating contour-plot of the difference between the analytic and numerical values at a couple of timesteps
for i in [0,int(len(t)/8),int(len(t)/4)]:
    filename = "analytisk_minus_numerisk,n=" + str(n) +",dt=" +str(dt) + ",dx=" +str(dx) + ",M=" + str(m_max) + ",T=" +str(t[i])
    plt.contourf(X,Y,u_all[i]-analytical_values[i])
    plt.colorbar()
    plt.savefig(dir+filename+".png")
    plt.clf()
