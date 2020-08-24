# About:

This code simulates a system of gravitating particles and plots their distribution in space (plots in 2D, although a 3D distribution of particles may be chosen). The updates to the positions and velocities of particles are calculated with the leapfrog method (similar to Euler’s method). 

In each update, the correlation function of the separation of pairs of particles and the corresponding power spectrum (the fourier transform of the correlation function) is computed. The correlation function is essentially a plot of the excess probability of finding a pair of particles separated by some distance, compared to a uniform distribution. The power spectrum is the same but in frequency space (it contains all the same information) - it tells you how likely you are to find a certain spatial frequency of particles. The particle distribution, correlation function and power spectrum are plotted and saved in a video to be viewed after the code has finished running.

## See this link to a youtube video explaining more about this code / demonstrating some interesting concepts I learn't while working on this project!

## Table of contents:
* [Running the code / Technologies required](#running-the-code---technologies-required)
* [Notes on the accuracy of the simulation and approximations used](#notes-on-the-accuracy-of-the-simulation,-approximations-used-and-limitations)
* [Equations/physics](#equations---physics)
* [Overview of hierarchy of functions](#overview-of-hierarchy-of-functions)
* [Interpreting the results](#tips-on-interpreting-the-results)

# Running the code / Technologies required
This code was created using Python 3.8. To run the python simulation code in your command prompt:
```
python Grav_particle_sim_AbbeW.py
```
Only the python file Grav_particle_sim_AbbeW.py is required. The modules numpy, matplotlib and scipy are required.

A description of all the input parameters required for the code is written in the python file - just call ```help(run_particle_sim())``` within the code to get a description of the input parameters printed to the terminal. Default values are set for each parameter, but the user has the option to enter them themselves.

Go to line 309 to change the file path the video is saved to and the name of the video:
```
ani.save(r"C:\Users\abbew\PythonProjects\Project1_PHYS4080\test.Mp4")
```
Note that you might see the following warning when you run the code:
```
RuntimeWarning: invalid value encountered in true_divide
  corr = np.where(unif_h != 0, ((500/Np)**2) *(h/unif_h) - 1, 0)  
```
You can ignore this warning, the function ```np.where()``` from the numpy module is used on this line to avoid dividing by zero (although it ironically spits out the warning), so this never actually occurs at all.

# Notes on the accuracy of the simulation, approximations used and limitations
The leapfrog method has been chosen over the Euler method because it is more accurate for the same amount of computational effort and is able to exactly serve the energy of the particles in the system and angular momentum. One should note the leapfrog/Euler method are not generally very efficient for very large numbers of particles. This is because the number of computations required to calculate the acceleration of each particle by each other partucke due to their gravitational interactions involves on the order of N(N-1) computations for each update of the particles velocities, where N refers to the number of particles in the simulation. The simulation will become slower as the number of particles is increased (it gets quite slow for more than a few hundred particles).
The particles require a ‘softening length’ to be included in the equation to calculate the acceleration by each other particle, because the particles are simulated as ‘points’ which is not representative of any real physical situations – point particles don’t exist in nature. If the softening length is not included to compensate for the fact the particles may reach a separation of 'zero', the acceleration of the particles by others can approach infinity (resulting in NaNs) if they get too close. The softening length serves to stop the acceleration becoming too big and suppresses close interactions. This approximation is at the sacrifice of accuracy and therefore the resolution of the simulation is limited by the softening length (no real physical processes can be considered accurate in the simulation on scales smaller than this length). The user is advised also, that for better accuracy, the timestep in each update should be chosen to be around Δt = ```final_t/Nt``` ~L/v  where L is the softening length, v is the typical velocity of particles in the system and ```final_t``` just refers to a variable that sets how long the simulation runs for.

Another points to consider is that the way in which the power spectrum is computed is not extremely accurate. The power spectrum is computed via a numerical integration of the correlation function, however the histogram for the correlation function that is plotted is generally not smooth (unless a very large number of particles and bins are used). The correlation function is interpolated first to try smooth it for the numerical integration to obtain the power spectrum, and gaussian quadrature integration is used to calculate the power spectrum, however the results are still not perfect and additionally the integral with limits from zero to infinity (see next section) is truncated since the space being simulated is only finite. The power spectrum output should be interpreted with this in mind.

A periodic boundary condition is encoded so that ‘new’ particles are always coming into the box when particles ‘exit’ the box. This is because only a finite length of space can be simulated. If the periodic boundary condition did not exist, particles may fly out of the box and none would be left for us to see. This boundary condition is equivalent to assuming the ‘universe’ being simulated is very large in its extent compared to the simulation volume and is approximately homogeneous and isotropic - so that the number of particles entering approximately equals those leaving the volume of space in consideration. 

# Equations/physics

Due to the finite length of the box, the correlation function and power spectrum are limited by a maximum possible distance particle pairs can be separated in space/minimum frequency and thus are truncated. The correlation function ξ(r) is calculated using the Peebles-Hauser approximation:

<img src="https://render.githubusercontent.com/render/math?math=\xi(r)=(\frac{N_r}{N_d})^2\frac{DD(r)}{RR(r)}-1"> 

RR(r) is the number of particles found in a uniform random distribution that are separated by some distance in an interval r+dr, DD(r) is the same for the particles in the simulation distribution, <img src="https://render.githubusercontent.com/render/math?math=N_r">  is the number of particles in the uniform distribution, <img src="https://render.githubusercontent.com/render/math?math=N_d ">  is the number of particles in the simulation. This function can easily be computed by simply counting the number of particles separated by various distances and binning them.

The power spectrum is computed via numerical integration of the correlation function. Note the integral in the simulation does not actually go to infinity because the finite length of the box truncates the correlation function at separations of 2√2L (L being half the box length). The power spectrum is given by:

<img src="https://render.githubusercontent.com/render/math?math=P(k)=2\pi \int_0^{\infty} x^2 \sinc{(kx)} dx  ">      

In general the force between two particles of mass i and j at separation r is computed as:

<img src="https://render.githubusercontent.com/render/math?math=F_{ij} = -\frac{G M_i M_j}{r}^{2}">

The softening length η (as discussed in the previous section) is incorporated into the equation below for the force between two particles of mass i and j at separation r when r is less than η:

<img src="https://render.githubusercontent.com/render/math?math=F_{ij} = -\frac{G M_i M_j}{r %2B \eta }^{2}">

These equations allow us to calculate the acceleration on each particle by all others.

The leapfrog method uses the following scheme to update the particles positions and velocities in each update (here i is refering to the ith update):
1)	Update the position of each particle x by half a timestep (only on the first update, v here is the velocity of the particle): 

<img src="https://render.githubusercontent.com/render/math?math=x_{\frac{1}{2}} =x_i %2B \frac{ \Delta t v }{2}">

2)	Calculating the update to the velocity of each particle by a full timestep (here a is the acceleration on the particle):

<img src="https://render.githubusercontent.com/render/math?math=v_{i%2B1} = v_i %2B \Delta t a_i">

3)	Update the position of each particle by a full timestep:
     
<img src="https://render.githubusercontent.com/render/math?math=x_{i%2B1%2B\frac{1}{2}} = x_{i%2B\frac{1}{2}} %2B \Delta t v_{i%2B1}">

Steps 2) to 3) are repeated in every update after the first one through the simulation.

# Overview of hierarchy of functions
The main program just gets user input and then calls ```run_particle_sim()```.

```run_particle_sim()```:
- set ups initial values of all variables based on input parameters to the function.
- sets up the settings and commands for plotting the particle distribution, correlation function and power spectrum.
- calls the function ```animation.FuncAnimation()``` which calls the function ```update()``` many times
to create each frame of the simulation video.

The function ```update()```:
- calculates updates to each particles' position and velocity (calls ```acceleration()``` to calculate the acceleration on each of the particles by all others and ```apply_boundary()``` to apply the boundary conditions after their positions have been updated) - it then sets the results computed for the particles positions' in the ```scatter()``` plot of the particle distribution via the command ```points```.
- calculates the correlation function by calling ```separation()``` which computes the separations between all the particles and bins them, then the Peeble-Hauser approximation is used to calculate the correlation function. Finally, it sets the results computed in the plot/histogram for the correlation function via the command ```points2```.
- calculates the power spectrum by calling ```power_spectrum2()``` (which calls ```Pk()```) to numerically integrate the correlation function to get the power spectrum using the equation for the fourier transform of the correlation function given previously. Then, it sets the results computed in the plot of the power spectrum via the command ```points3```.

# Tips on interpreting the results 
## The particle distribution

*	if the particles just fly out, there initial velocity or their masses may be very high, so they are accelerated away from each after initially maybe collapsing towards each other. ALthough they should be attracted to each other due to their gravity, they are able to gain a lot of acceleration and be flung out of the n-body system as they gain a lot of kinetic energy

* if they don’t move at all, they may be too far apart and or have an initial velocity that is too low (the gravitational attraction between the particles is too weak to make them move much over the time during the simulation) or perhaps they have very small masses

*	if they just fly around and don’t seem to interact, the velocity is probably very high (their kinetic energy is greater than the potential energy due to gravitational attractions) and thus there appears to be no evidence of any gravitational attraction

## Correlation function:
*	if the correlation function is very big at small separations (r), many particles should be close to each other in the plot showing the particle distribution (we expect the probability of finding particles at small separations to be large!) - the same idea in the opposite situation

* If the plot is totally zero everywhere the particles have no correlations at all - so the distribution in the particle distribution plot should quite uniform

## Power spectrum:
* If the power spectrum is large at low frequencies (k) the particles may have non-zero correlations at larger separations, and vice versa (you are likely to find particles separated by small distances in the particle distribution if the power spectrum is greater at large frequencies)

Overall, we expect that the results for the correlation function and power spectrum should have some reflection on the distribution of particles (although as previously explained, the power spectrum results need to be interpreted carefully due to the imperfect method being employed to determine the power spectrum). If the distribution ‘t’ (see descriptions of input parameters to the function and example video link) is input to the simulation, the correlation function will for example, initially have peaks at three main places/separations: at very small separations due to the four clumps of close-by particles, at separations equal to the length/side of the ‘box’ the 4 blobs of particles make and at another separation corresponding to the diagonal of the box the blobs of particles will make when using this initial setting of the particles distribution in space.
