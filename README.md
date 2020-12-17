# Project 5
I will denote Crank-Nicolson as CN, Backward-Euler as BE and Forward-Euler as FE in the Readme.

To run the programs, go into the c++ directory of the makefile and write "make all" in the terminal to compile and execute. Or "make compile", "make execute".

This requires that you have a folder "textfiles" in the main directory. In the makefile you have a line called conditions, this can be changed to run the simulation with different variables. The first argument is the size "n" of the vector u(x,t_j) at a spesific timestep j. The second is the numer of timesteps the simulations will iterate over. The third is the number of timesteps saved to file. 0 for saving all timesteps, 1 for saving every other, 2 for saving every third ect. The fourth argument is alpha. The fifth is whether to run the simulation in one or two spacial dimentions. "1d" for one dimention and "2d" for two dimentions.

This will only create the textfiles. The files in the folder "python" will is used to analyse the data in the textfiles. None of the python files takes any input arguments, so they can be run with the command "python3 filename", where "filename" is the name of the python file. To run the files you need to create a folder "figures" in the main directory.

The file "2d_plot.py" will save two figures for two different timesteps with BE, FE and CN plotted on top of eachother. The timestep numers are on line 34 and can be changed to whatever timesteps you want to look at (1,2,3,...) timesteps can be chosen. 

The file "3d_plot.py" will save one figure for BE, one for FE and one for CN. The plots will be a 3d-plot of x,u(x) and the time t, so we can see how the different models behave in time.

The file "analytic.py" sets up an analytic solution for u(x,t), and saves a plot of the analytic solution against BE, FE and CN.

The file "analytic_mot_crank.py" is used to check different "m_max", where "m_max" is the maximum value that we let the sum in the analytic solution go to, since we to approximate the analytic solution. This is then plotted against the CN model. These values can be changed on line 42.

The file "error_func.py" saves a figure of the error epsilon(t) of the numerical solution against the analytic solution for BE,FE,CN.

The file "analytic_2dim.py" is for the 2-dimensionional case. This will first calculate analytic values. Then save a plot of the error epsilon(t). Then save 3d plots for a number of given timesteps for the analytic and numerical solution. These timesteps are chosen at line 117 and 121. Last it will save contour plots of the numerical minus the analytical solution for a number of given timesteps. In other words the difference between the analytic and numerical solution. These timesteps are chosen at line 129.
