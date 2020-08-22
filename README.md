# About:

This code simulates a system of gravitating particles and plots their distribution in space (plots in 2D, although a 3D distribution may be chosen). The updates to the positions and velocities of particles are calculated with the leapfrog method (similar to Euler’s method). 

In each update, the correlation function of the separation of pairs of particles and the corresponding power spectrum (the fourier transform of the correlation function) is computed. The correlation function is essentially a plot of the excess probability of finding a pair of particles separated by some distance, compared to a uniform distribution. The power spectrum is the same but in frequency space (it contains all the same information) - it tells you how likely you are to find a certain spatial frequency of particles. The particle distribution, correlation function and power spectrum are plotted and saved in a video to be viewed after the code has finished running.

# Example video:



# Running the code

Only the python file Grav_particle_sim_AbbeW.py is required. The modules numpy, matplotlib and scipy are required.

A description of all the input parameters to the main function in the code ( run_particle_sim() ) is written in the python file (just call help(run_particle_sim()) within the code to get a description of the input parameters printed to the terminal, scroll to bottom of code and uncomment line 317). Default values are set for each parameter, but the user has the option to enter them themselves.

Go to line 312 to change the file path the video is saved to and the name of the video:
```
ani.save(r"C:\Users\abbew\PythonProjects\Project1_PHYS4080\test.Mp4")
```
You might see the following warning when you run the code:
```
RuntimeWarning: invalid value encountered in true_divide
  corr = np.where(unif_h != 0, ((500/Np)**2) *(h/unif_h) - 1, 0)  
```
You can ignore this warning, the function np.where() is used on this line to avoid dividing by zero (although it ironically spits out the warning), so this never actually occurs at all.

# Accuracy of the simulation and approximations used

# Equations/physics

# Brief explanation of function hierarchy

# Interpreting the results 
