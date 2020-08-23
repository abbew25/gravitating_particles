# About:

This code simulates a system of gravitating particles and plots their distribution in space (plots in 2D, although a 3D distribution of particles may be chosen). The updates to the positions and velocities of particles are calculated with the leapfrog method (similar to Euler’s method). 

In each update, the correlation function of the separation of pairs of particles and the corresponding power spectrum (the fourier transform of the correlation function) is computed. The correlation function is essentially a plot of the excess probability of finding a pair of particles separated by some distance, compared to a uniform distribution. The power spectrum is the same but in frequency space (it contains all the same information) - it tells you how likely you are to find a certain spatial frequency of particles. The particle distribution, correlation function and power spectrum are plotted and saved in a video to be viewed after the code has finished running. 

## Table of contents:
* [Running the code / Technologies required](# Running the code / Technologies required)
* [Notes on the accuracy of the simulation and approximations used](# Notes on the accuracy of the simulation and approximations used)
* [Equations/physics](# Equations/physics)
* [Overview of hierarchy of functions](# Overview of hierarchy of functions)
* [Interpreting the results](# Interpreting the results)

# Running the code / Technologies required
This code was created using Python 3.8. To run the python simulation code in your command prompt:
```
python Grav_particle_sim_AbbeW.py
```
Only the python file Grav_particle_sim_AbbeW.py is required. The modules numpy, matplotlib and scipy are required.

A description of all the input parameters required for the code is written in the python file - just call help(run_particle_sim()) within the code to get a description of the input parameters printed to the terminal. Default values are set for each parameter, but the user has the option to enter them themselves.

Go to line 312 to change the file path the video is saved to and the name of the video:
```
ani.save(r"C:\Users\abbew\PythonProjects\Project1_PHYS4080\test.Mp4")
```
Note that you might see the following warning when you run the code:
```
RuntimeWarning: invalid value encountered in true_divide
  corr = np.where(unif_h != 0, ((500/Np)**2) *(h/unif_h) - 1, 0)  
```
You can ignore this warning, the function np.where() is used on this line to avoid dividing by zero (although it ironically spits out the warning), so this never actually occurs at all.

# Notes on the accuracy of the simulation and approximations used
The leapfrog method has been chosen over the Euler method because it is more accurate for the same amount of computational effort and is able to exactly serve the energy of the particles in the system and angular momentum. One should note the leapfrog/Euler method are not efficient for large numbers of particles as the number of computations to calculate the acceleration of each particle by each other due to the gravitational interactions as it is a direct n-body method that involves N(N-1) computations for each update of the particles velocities, where N is the number of particles in the simulation. The simulation will become slower as the number of particles is increased (it gets quite slow for more than a few hundred particles).
The particles require a ‘softening length’ to be included in the equation to calculate the acceleration by each other particle, because the particles are simulated as ‘points’ which is not representative of any real physics – point particles don’t really exist in nature. If the softening length is not included to compensate for the fact the particles may reach a separation of 'zero', the acceleration of the particles by others can approach infinity (resulting in NaNs) if they get too close. The softening length stops the acceleration becoming too big and suppresses close interactions. This approximation is at the sacrifice of accuracy and therefore the resolution of the simulation is limited by the softening length (no real physical processes can be considered accurate in the simulation on scales smaller than this length). The user is advised also that for accuracy the timestep in each update should be chosen to be around Δt = final_t/Nt ~L/v  where L is the softening length and v is the typical velocity of particles in the system.

Another points is that the way in which the power spectrum is computed is not extremely accurate. The power spectrum is computed via a numerical integration of the correlation function using the trapezoid method, however the histogram for the correlation function that is plotted is generally not smooth (unless a very large number of particles and bins are used). The correlation function is interpolated to try smooth it for the numerical integration of the power spectrum. However you can see in the example videos I have linked that the power spectrum still looks very not smooth! This is a result of integration of the function shown in the next section.

A periodic boundary condition is encoded so that ‘new’ particles are always coming into the box when particles ‘exit’ the box. This is because only a finite length of space can be simulated. If the boundary condition did not exist, particles may fly out of the box and none would be left to see. This boundary condition is equivalent to assuming the ‘universe’ being simulated is very large in its extent compared to the simulation volume and is approximately homogeneous and isotropic. 

# Equations/physics

Due to the finite length of the box, the correlation function and power spectrum are limited by a maximum possible distance particle pairs can be separated /minimum frequency in space and thus are truncated. The correlation function ξ(r) is calculated using the Peebles-Hauser approximation:

<img src="https://render.githubusercontent.com/render/math?math=\xi(r)=(\frac{N_r}{N_d})^2\frac{DD(r)}{RR(r)}-1"> 

RR(r) is the number of particles found in a uniform random distribution which are separation by some distance in the interval r+dr, DD(r) is the same for the particles in the simulation distribution, N_r is the number of particles in the uniform distribution, N_d is the number of particles in the simulation.

The power spectrum is computed via numerical integration. Note the integral in the simulation does not actually go to infinity because the finite length of the box truncates it at 2√2 L (L being half the box length). The power spectrum is given by:

<img src="https://render.githubusercontent.com/render/math?math=P(k)=2\pi \int_0^{\infty} x^2 \sinc{(kx)} dx  ">      

The softening length η is incorporated into the equation below for the force  between two particles of mass i and j at separation r when r is less than η:

<img src="https://render.githubusercontent.com/render/math?math=F_{ij} = -\frac{G M_i M_j}{r %2B \eta }^{2}">

The leapfrog method uses the follow scheme to update particle positions and velocities does the following in each update (i refering to the ith update):
1)	Update the position of particles by half a timestep (only on the first update): 

<img src="https://render.githubusercontent.com/render/math?math=x_{\frac{1}{2}} =x_i %2B \frac{ \Delta t }{2v}">

2)	Calculating the update to the velocity of the particles by a full timestep:

<img src="https://render.githubusercontent.com/render/math?math=v_{i%2B1} = v_i %2B \Delta t a_i">

3)	Update the position of the particle by a full timestep:
     
<img src="https://render.githubusercontent.com/render/math?math=x_{i%2B1%2B\frac{1}{2}} = x_{i%2B\frac{1}{2}} %2B \Delta t v_{i%2B1}">

Steps 2) to 3) are repeated in every update after the first one.

# Overview of hierarchy of functions
The main program just gets user input and then calls ```run_particle_sim()```.

```run_particle_sim()```:
- set ups initial values of parameters /quantities needed / initialises variables
- sets up the settings and commands for plotting
- calls the function ```animation.FuncAnimation()``` which calls the function ```update()``` many times

The function ```update()```:
- calculates updates to each of the particles in the distribution (calls ```acceleration()``` to calculate the acceleration on each of the particles and ```apply_boundary()``` to apply the boundary condition). It then sets the data in the ```scatter()``` plot for the particle distribution for the update with the command ```points``` 
- calculates the correlation function by calling ```separation()``` and then sets the data in the plot for that update with the command ```points2```
- calculates the power spectrum by calling ```power_spectrum2()``` (which calls ```Pk()```) to numerically integration the correlation function to get the power spectrum and then sets the data in the plot for that update with the command ```points3```

# Interpreting the results 
1) The particle distribution:

•	if the particles just fly out, there initial velocity or their masses may be very high, so they are accelerated away from each. ALthough they should be attracted to each other due to their gravity, they are able to gain a lot of acceleration and be flung out of the n-body system as they gain a lot of energy

•	if they don’t move at all, they may be too far apart and or have an initial velocity that is too low (the gravitational attraction between the particles is too weak to make them move much over the time during the simulation)

•	if they fly around and don’t seem to interact, the velocity is probably very high (their kinetic energy is greater than the potential energy due to gravitational attractions)

2) Correlation function:
•	if the distribution is very big at small r, many particles should be close to each other in the plot showing the particle distribution (we expect the probability of finding particles at small separations to be large!) The same idea in the opposite situation

•	If the plot is totally zero everywhere the particles have no correlations at all - so the distribution in the particle distribution plot should appear uniform!

3) Power spectrum:
•	If the power spectrum is large at low frequencies (k) the particles may have non-zero correlations at larger separations, and vice versa (you are likely to find particles separated by small frequencies or large etc.)

Overall, the results for the correlation function and power spectrum should have some reflection on the distribution of particles. If the distribution ‘t’ is input to the simulation, the correlation function will for example, initially have peaks like a delta function at 3 places: at very small separations due to the four clumps of close by particles, at separations equal to the length/side of the ‘box’ the 4 blobs of particles make and at another separation corresponding to the diagonal of the box the blobs of particles will make.
