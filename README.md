# About:

This code simulates a system of gravitating particles and plots their distribution in space (plots in 2D, although a 3D distribution of particles may be chosen). The updates to the positions and velocities of particles are calculated with the leapfrog method (similar to Euler’s method). 

In each update, the correlation function of the separation of pairs of particles and the corresponding power spectrum (the fourier transform of the correlation function) is computed. The correlation function is essentially a plot of the excess probability of finding a pair of particles separated by some distance, compared to a uniform distribution. The power spectrum is the same but in frequency space (it contains all the same information) - it tells you how likely you are to find a certain spatial frequency of particles. The particle distribution, correlation function and power spectrum are plotted and saved in a video to be viewed after the code has finished running.

## You can follow this link to a youtube video related to this code, which very briefly demonstrates some interesting concepts I learn't while working on this project.
<iframe width="560" height="315" src="https://www.youtube.com/embed/9h1CMGiAvKk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

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

A description of all the input parameters required for the code is written in the python file - just call ```help(run_particle_sim())``` within the code to get a description of the input parameters printed to the terminal. Default values are set for each parameter so you can just run the code quickly, but the user has the option to enter them themselves.

You might like to go to line 309 to change the file path the video is saved to and the name of the video:
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
The leapfrog method has been chosen over the Euler method because it is more accurate for the same amount of computational effort and is able to exactly serve the energy of the particles in the system and angular momentum. One should note however, that the leapfrog/Euler method are not generally very efficient for very large numbers of particles. This is because the number of computations required to calculate the acceleration of each particle by each other partucke due to their gravitational interactions involves on the order of many computations for each update of the particles velocities. This simulation will become slower quite quickly as the number of particles is increased (it gets quite slow for more than a few hundred particles). Therefore it's performance is limited to simulating just a few hundred particles.
The particles require a ‘softening length’ as an input parameter to the code so it can be included in the equations used to calculate the acceleration on every particle by each other. This is because the particles are simulated as ‘points’ which is not representative of any real physical situations – point particles don’t exist in nature. If the softening length is not included to compensate for the fact the particles may reach a separation of 'zero', the acceleration of the particles by each other may approach infinity (resulting in NaNs) if they get too close. The softening length serves to stop the acceleration becoming too big and suppresses close interactions between the particles. This approximation is at the sacrifice of accuracy on smaller scales and therefore the resolution of the simulation is limited by the softening length chosen (no real physical processes can be considered accurate in the simulation on scales smaller than this length). The user is advised also, that for the best accuracy, the timestep in each update should be chosen to be around Δt = ```final_t/Nt``` ~L/v  where L is the softening length, v is the typical velocity of particles in the system. ```final_t``` just refers to a variable that sets how long the simulation runs for in Gyrs, while ```Nt``` is the number of updates that occur, thus together they set Δt.

Another important point to consider is that the way in which the power spectrum is computed is unfortunately not extremely accurate. The power spectrum is computed via a numerical integration of the correlation function via a fourier transform, however the correlation function that is calculated from the particle distribution is generally not smooth (unless a very large number of particles and bins are used, but even then it isn't very smooth). The correlation function is interpolated first to try smooth it for the numerical integration to obtain the power spectrum, and fixed gaussian quadrature integration is used to calculate the power spectrum. However, the results are still found not be extremely accurate. Additionally, the integral over the correlation function has limits from zero to infinity (see next section for more detail) but is has to be truncated since the space being simulated is only finite. The power spectrum output should be interpreted with this in mind.

A periodic boundary condition is encoded so that ‘new’ particles are always coming into the box when particles ‘exit’ the box. This is because only a finite length of space can be simulated. If the periodic boundary condition did not exist, particles may fly out of the box and none would be left for us to see. This boundary condition is equivalent to assuming the ‘universe’ being simulated is very large in its extent compared to the simulation volume and is approximately homogeneous and isotropic - so that the number of particles entering approximately equals those leaving the volume of space in consideration. 

# Equations/physics

Due to the finite length of the box, the correlation function and power spectrum are limited by a maximum possible distance particle pairs can be separated in space/minimum frequency and thus are truncated. The correlation function ξ(r) is calculated using the Peebles-Hauser approximation:

<img src="https://render.githubusercontent.com/render/math?math=\xi(r)=(\frac{N_r}{N_d})^2\frac{DD(r)}{RR(r)}-1"> 

RR(r) is the number of particles found in a uniform random distribution that are separated by some distance in an interval r+dr, DD(r) is the same for the particles in the simulation distribution, <img src="https://render.githubusercontent.com/render/math?math=N_r">  is the number of particles in the uniform distribution, <img src="https://render.githubusercontent.com/render/math?math=N_d ">  is the number of particles in the simulation. This function can easily be computed by simply counting the number of particles separated by various distances and binning them to get RR(r) for the particle distribution and DD(r) for a uniform distribution.

The power spectrum is computed via numerical integration of the correlation function. Note the integral in the simulation does not actually go to infinity because the finite length of the box truncates the correlation function at separations of 2√2L (L being half the box length). However, the power spectrum is given by:

<img src="https://render.githubusercontent.com/render/math?math=P(k)=2\pi \int_0^{\infty} x^2 \sinc{(kx)} dx  ">   

In the code, the integral is taken from just zero to the maximum separation any particles can have in the box, and for larger separations the correlation function is assumed to be zero and thus has no contribution to the power spectrum. This means we are assuming that there are no correlations on larger length scales than the box. For the integral the correlation function is interpolated by using the scipy function splrep and passing the spline is returns to integral above to compute the correlation at any given x, allowing the integration to be smoother. Then the result of integration is interpolated again so a smooth power spectrum can be plotted.

In general the force between two particles of mass i and j at separation r is computed as:

<img src="https://render.githubusercontent.com/render/math?math=F_{ij} = -\frac{G M_i M_j}{r}^{2}">

This code calculates the force on each individual particle by every single other particle via the above equation above one at a time. This can be done because the positions of all particles in every dimension is stored in a vector that is updated in each step of the simulation. 
The softening length η (as discussed in the previous section) is incorporated into the equation below for the force between two particles of mass i and j at separation r when r is less than η:

<img src="https://render.githubusercontent.com/render/math?math=F_{ij} = -\frac{G M_i M_j}{r %2B \eta }^{2}">

Thus even if r = 0 for some two particles, the force between the two particles will not be infinite as long as η is not chosen to be zero. 

The above equations allow us to calculate the acceleration on each particle by all others.

The leapfrog method uses the following scheme to update the particles positions and velocities in each update (in the equations below, i is refering to the ith update):
1)	Update the position of each particle x by half a timestep (only on the first update, v here is the velocity of the particle): 

<img src="https://render.githubusercontent.com/render/math?math=x_{\frac{1}{2}} =x_i %2B \frac{ \Delta t v }{2}">

2)	Calculating the update to the velocity of each particle by a full timestep (here a is the computed acceleration on the particle by all other particles on this step):

<img src="https://render.githubusercontent.com/render/math?math=v_{i%2B1} = v_i %2B \Delta t a_i">

3)	Update the position of each particle by a full timestep (here v is the newly updated velocity from step 2)):
     
<img src="https://render.githubusercontent.com/render/math?math=x_{i%2B1%2B\frac{1}{2}} = x_{i%2B\frac{1}{2}} %2B \Delta t v_{i%2B1}">

Steps 2) to 3) are repeated in every update after the first one over the simulation so the motion of the particles can be seen.

# Overview of hierarchy of functions
The main program just gets user input and then calls ```run_particle_sim()```.

```run_particle_sim()```:
- set ups initial values of all variables based on input parameters to the function.
- sets up the settings and commands for plotting the particle distribution, correlation function and power spectrum.
- calls the function ```animation.FuncAnimation()``` which calls the function ```update()``` many times
to create each frame of the simulation video using various functions.

The function ```update()```:
- calculates updates to each particles' position and velocity (calls ```acceleration()``` to calculate the acceleration on each of the particles by all others and ```apply_boundary()``` to apply the boundary conditions after their positions have been updated) - it then sets the results computed for the particles positions' in the ```scatter()``` plot of the particle distribution via the command ```points```.
- calculates the correlation function by calling ```separation()``` which computes the separations between all the particles and bins there separations, then the Peeble-Hauser approximation is used to calculate the correlation function. Then it sets the results computed in a plot/histogram for the correlation function via the command ```points2```.
- it finally calculates the power spectrum by calling ```power_spectrum2()``` (which calls ```Pk()```) to numerically integrate the correlation function (computed in the previous step) to get the power spectrum using the equation described previously for the fourier transform of the correlation function. Then, it sets the results computed in the plot of the power spectrum via the command ```points3```.

# Tips on interpreting the results 
## The particle distribution

*	if the particles just fly out or away from each other, there initial velocity or their masses may be very high, so they are accelerated away from each after initially maybe collapsing towards each other. ALthough they should be attracted to each other due to their gravity, they are able to gain a lot of acceleration and be flung out of the n-body system as they gain a lot of kinetic energy 

* if they don’t move at all, they may be too far apart and or have an initial velocity that is too low (the gravitational attraction between the particles is too weak to make them move much over the time during the simulation) or perhaps they have very small masses and thus there gravitational attraction is too weak

*	if they just fly around and don’t seem to interact, the velocity is probably very high (their kinetic energy is greater than the potential energy due to gravitational attractions) and thus there appears to be no evidence of any gravitational attraction

* by playing around with the input parameters for the size of the box, the initial particle distribution, the masses of the particles and initial velocities, it is possible to find a situation in which the particles may collapse into little clusters and gravitate

## Correlation function:
*	if the correlation function is very big at small separations (r), many particles should be close to each other in the plot showing the particle distribution (we expect the probability of finding particles at small separations to be large!) - the same idea applies in the opposite situation

* If the plot is totally zero everywhere the particles may have no correlations at all - so the distribution in the particle distribution plot should appear quite uniform, because this means there distribution is totally random

## Power spectrum:
* If the power spectrum is large at low frequencies (k) the particles may have non-zero correlations at larger separations, and vice versa (you are likely to find particles separated by small distances in the particle distribution if the power spectrum is greater at large frequencies)

* overall the power spectrum contains the same information as the correlation function but only differs in that it is showing the correlations between particle positions in frequency space

Overall, we expect that the results for the correlation function and power spectrum should have some reflection on the distribution of particles (although as previously explained, the power spectrum results need to be interpreted carefully due to the imperfect method being employed to determine the power spectrum in this code). If the distribution ‘t’ (see descriptions of input parameters to the function) is input to the simulation, the correlation function will for example, initially have peaks at three main places/separations: at very small separations due to the four clumps of close-by particles, at separations equal to the length/side of the ‘box’ the 4 blobs of particles make and at another separation corresponding to the diagonal of the box the blobs of particles will make when using this initial setting of the particles distribution in space.
