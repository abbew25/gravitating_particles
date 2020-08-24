# Abbe Whitford PHYS4080

# About:
# This code simulates gravitating particles and also calculates their correlation
# function and power spectrum for the particle distribution for each time step.
# Change line 309 to in ani.save() before you run this code to save the video it generates
# with plots in a folder you would like it to go in.

# Input files: None. 

# Output files: Save mp4 file called 'test.mp4' of the particle simulation (video includes
# plot of distribution of particles in space, their correlation function and power spectrum)

# Scroll to the end of the code to change the default input parameters to run_particle_sim().
# The program also asks for user input - you may choose to specify parameters each time you
# run the code instead.

# Imports: -----------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as m
from scipy import interpolate, integrate
#---------------------------------------------------------------------------------------------

def run_particle_sim(particle_mass_max = 10, frame_duration = 50, t_final = 1000, 
L = 400, num_bins = 500, softening_length = 20, Np = 30, Nt = 50, Nd = 2, v_max = 0.010, 
n = 30, dist_type = 'g'):
    """4
    Input parameters are (in order):
        particle_mass_max - maximum mass set for any particle (in solar mass units)

        frame_duration - how long each frame runs for in video (just speeds up or slows down 
        video)

        t_final - the time the simulation ends (myrs) at (from 0 to t_end = e.g. 1000 Myrs)

        L - half of the height and width of the 2d box scatterplot

        num_bins - the number of bins used in the correlation function and power spectrum 
        (sets smoothing scale)

        softening_length - the softening length (sets resolution scale)

        Np - number of particles (number of updates /frames computed)

        Nt - number of time steps/updates (sets deltat)

        Nd - number of spatial dimensions (pick 2 or 3 only)

        v_max - max. initial velocity of any particle (solar masses/Myrs)

        n - n is the standard deviation of a gaussian distribution (if this OR the blobs
        -of-particles type distributions are chosen for the initial positions of particles)

        dist_type  - the type of distribution out of uniform 'u', gaussian 'g' or funky
        blobs 'f' or 't' ('f' is many spots with blobs of particles in guassian distributions 
        along a diagonal, while 't' is 4 blobs of particles in a square, each blob has a gaussian 
        distribution also). Pick a number of particles divisible by 10 if using 'f' and a number 
        of particles divisible by 4 if using 't'.
    """


    #---------------------------------------------------------------------------------------------
    time = 0                          # time sim 'starts'
    # give particles random masses up to some max mass.
    particle_masses = particle_mass_max - 0.5*particle_mass_max*np.random.random((1,Np))

    # give particles different transparency of color in scatterplot according to mass
    # (using same vector also to change marker size of particles in plot according to mass)
    colors_by_mass = np.vstack((       # red, green, then blue amounts of color for each particle 
    (particle_masses[0,:] ),
    (particle_masses[0,:]*0.5 ), 
    (particle_masses[0,:]*0.1 ))).T 
            
    # reshape and normalize colors_by_mass to pass to scatterplot()
    colors_by_mass = colors_by_mass.reshape((Np,3))/particle_mass_max
    np.random.seed(3000)               # For reproducibility, set a seed for  
                                       # randomly generated inputs. 
    delta_t = (t_final)/(Nt)           # sets amount of time passed in each time frame # (Myrs)
    G = 0.0044                         # gravitational constant in units of 
                                       # pc^3 /solar mass * Myrs^2

    # maximum seperation considered - setting up for correlation function histogram
    max_sep = 2*np.sqrt(2)*L           # just initializing variable            
    if Nd == 3:
        max_sep = 2*m.sqrt(3)*L
    smoothing_length = int(2*np.sqrt(2)*L/num_bins) # width of bins in histogram (correlation function)
    sbins = np.linspace(0,max_sep, num_bins)        # define bins for correlation function
    sbins[0] = 1e-6                                 # don't want frequency at infinity later.
    #set up stuff for frequency space histogram (power spectrum)
    k_max = 2*np.pi/(0.5*softening_length) # bit smaller than smallest length scale info. can be resolved on 
    k_min = 2*np.pi/(max_sep)          # largest length scale in simulation
    # Set initial velocities to be random fractions of the maximum (pc/Myrs)
    velocity = v_max*(1-2*np.random.random((Nd,Np)))
    position = np.zeros((Nd, Np))      # initialise array of positions of particles
    
    # set initial distribution of particles based on input to function
    if dist_type == 'u':
        # initial uniform distribution
        position = L - 2*L*np.random.random((Nd,Np))
    elif dist_type == 'g':   
        #guassian initial distribution
        position = np.random.normal(0, n, (Nd,Np))
    elif dist_type == 'f':
        #Set initial positions at spotty blobs within box along the diagonal
        p = int(Np/10)                     # number particles at each blob (for 'f' distribution)
        position = np.hstack( (np.random.normal(-750, n, (Nd,p)), 
                               np.random.normal(750, n,  (Nd,p)),
                               np.random.normal(500, n,  (Nd,p)), 
                               np.random.normal(-500, n, (Nd,p)), 
                               np.random.normal(250, n,  (Nd,p)), 
                               np.random.normal(-250, n, (Nd,p)), 
                               np.random.normal(-100, n, (Nd,p)), 
                               np.random.normal(-100, n, (Nd,p)), 
                               np.random.normal(0, n,  (Nd,p)), 
                               np.random.normal(0, n,  (Nd,p)) 
                            ) )
    elif dist_type == 't':
        # set 4 blobs in a square shape 
        p = int(Np/4)                     # number particles at each blob (for 't' distribution)
        g1 = np.random.normal(0, n, (Nd, p))
        g2 = np.random.normal(0, n, (Nd, p))
        g3 = np.random.normal(0, n, (Nd, p))
        g4 = np.random.normal(0, n, (Nd, p))
        g1[0] += L/2
        g2[0] -= L/2
        g1[1] += L/2
        g2[1] += L/2
        g3[0] -= L/2
        g4[0] += L/2
        g3[1] -= L/2
        g4[1] -= L/2
        position = np.hstack( ( 
                    g1, g2, g3, g4
                                 ) )
    else:
        position = np.random.normal(0, n, (Nd,Np))

    # end of setting up variables------------------------------------------------------
        

    #create/set up figure for simulation-----------------------------------------------
    # Set the axes on which the points will be shown
    plt.ion()                                           # Set interactive mode on
    fig = plt.figure(figsize=(16,8))     
    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95,wspace=0.15,hspace=0.3)
    
    # plotting gravitating particles        
    ax1 = plt.subplot(121, facecolor = 'black')         # Set the axes as the only set 
    ax1.set_xlim(-L,L)                                  # Set x-axis limits
    ax1.set_ylim(-L,L)                                  # Set y-axis limits
    ax1.set_xlabel('parsecs (pc)')
    ax1.set_ylabel('parsecs (pc)')
    ax1.set_title(r"Gravitating  particles (%d) with $m_{max.}$ = %d $m_{\odot}$." % (Np, 
    particle_mass_max))
    ax1.text(-L+L/5, L-L/5, "time = %.2f Myrs" % (time), c = 'white')
    
    
    # plotting histogram (correlation function)
    ax2 = plt.subplot(222)
    ax2.set_xlabel('Seperation r (pc) (Smoothing length = %.2f pc, max. box sep. = %.2f pc)' 
    % (smoothing_length, max_sep))
    ax2.set_ylabel('Correlation of pairs seperated by r (pc)')
    ax2.set_title("Correlation Function")
    ax2.set_xlim([softening_length, max_sep])
    ax2.set_ylim([0.01, 1000])
    plt.yscale("log")
    
    #plotting power spectrum 
    ax3 = plt.subplot(224)
    ax3.set_xlabel('k $pc^{-1}$')
    ax3.set_ylabel('P(k) ')
    ax3.set_title("Power Spectrum")
    ax3.set_xlim([k_min, k_max])
    ax3.set_ylim([0.0001, 10e10])
    plt.xscale("log")
    plt.yscale("log")
    
    # Create commands which will plot the positions of the particles distribution, correlation 
    # function and power spectrum
    # (commands are points, points2, points3 in that order)
    new_position = np.reshape(np.concatenate((position[0,:], position[1,:])), (2,Np)).T
    points = ax1.scatter(new_position[:,0], new_position[:,1], c = colors_by_mass, 
    s = 10*colors_by_mass[0,:])
    points2, = ax2.plot([],[],drawstyle='steps-post', c = 'red')
    points3, = ax3.plot([],[],drawstyle='steps-post', c = 'orange')

    pos_Random = L - 2*L*np.random.random((Nd, 500))    # get random distribution of particles 
    # (above variable to be used later in 'update' function)
    
    # define nested functions:--------------------------------------------------------
    
    
    # Create a function to apply periodic boundary conditions (just requires positions p of particles)
    def apply_boundary(p):
      newp = p
      newp[newp > L] = newp[newp > L]-2*L
      newp[newp < -L] = 2*L+newp[newp < -L]
      return newp  

    # function to calculate particle accelerations (just requires positions p of particles)
    def acceleration(p):                                            
        acc = np.zeros((Nd,Np))                                     #acceleration matrix to be returned
        #loop through the ith particle to get its acceleration
        for i in np.arange(Np): 
            difpos = p[:,:] - p[:,None,i]                           # distance between ith particles and all others
            copy_difpos = difpos.tolist()                           # need to remove ith particle (distance of ith particle from itself) from list
            for j in np.arange(Nd):                                 # converting back to list to make removal easy then convert back to numpy array
                copy_difpos[j].pop(i)
            difpos = np.array(copy_difpos) 
            abs_rs = np.zeros(Np-1).T                               # initialise and calculate |r_ij..| = sqrt( x^2 + y^2 ..) for each particle etc. 
            for j in np.arange(Nd):
                abs_rs += np.array(difpos[j,:]**2)
            abs_rs = np.sqrt(abs_rs)                                
            abs_rs[abs_rs < softening_length] += softening_length   # correct |r| with softening length if needed
            needed_particle_masses = np.delete(particle_masses, i)  # get rid of ith particle mass from list of masses
            x = (difpos*needed_particle_masses.T)*G/(abs_rs**3).T   # calculate the acceleration by each particle using their individual masses and 
                                                                    # seperations
            acc[:,i] +=  np.sum(x, axis=1).T                        
            
        return acc
    
   
    # Define function for procedure to update positions at each timestep (i is ith update)
    def update(i):
        # 1) Get positions and velocities and time variables
        nonlocal position,velocity                          
        nonlocal time          
        time += delta_t  # update time in this step

        # 2) Using the 'leapfrog' method (similar to Euler method) to update positions and velocities 
        # (updates are out of sync by half a timestep compared to Euler method):
        if (i == 0):
            position += velocity*delta_t/2                  # Increment positions with velocities
            position = apply_boundary(position)             # Apply boundary condition  
        velocity += acceleration(position)*delta_t        # increment velocity using acceleration
        position += velocity*delta_t                        # Increment positions with updated velocities                              
        position = apply_boundary(position)                 # Apply boundary condition


        # 3) Show 2D projection of first 2 position coordinates
        points.set_offsets(np.reshape(np.concatenate((position[0,:], position[1,:])), (2,Np)).T) 
        for txt in ax1.texts:
            txt.set_visible(False) # getting rid of 'time = something' text put on plot from previous step
        ax1.text(-L+L/5, L-L/5, "time = %.2f Myrs" % (time), c = 'white') #update current time in simulation on the plot
        

        # 4) get the correlation function and plot it
        h,x = np.histogram(np.ravel(np.tril(separation(position))),bins=sbins)  # get histogram for particle seperations
        unif_h, x = np.histogram(np.ravel(np.tril(separation(pos_Random))),bins=sbins) # get histogram for particle seperations (uniform dist.)
        # pos_Random was defined earlier in the code (to avoid generating random distribution repeatedly)
        # peebles-hauser approximation of correlation function:
        corr = np.where(unif_h != 0, ((500/Np)**2)*(h/unif_h) - 1, 0)  # calculate correlation function
        points2.set_data(x[:-1],corr) 
        

        # 4) get the power spectrum and plot it
        ks, Pks = power_spectrum2(x[:-1], corr)
        points3.set_data(ks,Pks) 
        
        return points, points2, points3
    
    # Function to find separations from position vectors (get DD(r) or RR(r)) (just requires positions of particles p)
    def separation(p):                  
        s = p[:,None,:] - p[:,:,None]   # find N x N x Nd matrix of particle seps.
        return np.sum(s**2,axis=0)**0.5 # return N x N matrix of scalar separations

    #function to integrate to get P(k) from numerical integration
    def Pk(x, k, tck):        # x value, k value, correlation(x) value       
        #pk = 2*np.pi*pow(x,2)*corr*(np.sinc(k*x))  
        pk = 2*np.pi*pow(x,2)*(interpolate.splev(x, tck, der=0))*(np.sinc(k*x)) 
        return pk

    # function to get matter power spectrum (just requires bin values (x-axis)
    #  and correlations (y - axis) from correlation function histogram)
    def power_spectrum2(bins, corrs):   
        # 1) get an interpolation of the correlation function to use for integration
        tck = interpolate.splrep(bins, corrs, k = 3, s=0) 

        # 2) Convert rs to ks (get separation bins in frequency space for integration and plotting)
        ks = 2*np.pi/bins
        ks = ks[::-1]           # get order of ks in increasing k

        # 3) Now use a for loop to get the P(k) (iterating through ks) 
        # by doing numerical integration with scipy
        # of the fourier transform of the correlation function
        P_ks = np.zeros(len(ks))
        for i in np.arange(len(ks)):    # iterate through ks to get P(k)
            # integrate for fixed k using gaussian quadrature
            integral = integrate.fixed_quad(Pk, 0, max_sep, (ks[i], tck), n = 400)
            P_ks[i] = integral[0]          # save result of integration to array

        # 4) get interpolation of the result of integration to try smooth the power spectrum
        tck = interpolate.splrep(ks, P_ks, k = 3, s=0) 
        ks = np.logspace(np.log10(k_min), np.log10(k_max), num = num_bins)
        P_ks = interpolate.splev(ks, tck, der=0)

        return ks, P_ks

    
    # end of nested functions:-------------------------------------------------
    
    # Create animation by calling update function!
    # https://matplotlib.org/api/_as_gen/matplotlib.animation.FuncAnimation.html
    ani = animation.FuncAnimation(fig, update, frames=Nt,interval = frame_duration)
    
    # save to mpeg file: 
    ani.save(r"C:\Users\abbew\PythonProjects\test3.Mp4")
    print('Finished.')
    return 0


#help(run_particle_sim())

# get user input to specify parameters or just run with simulation with default parameters:
res = input("Would you like to use the default parameters set for this code? Enter yes or no: ")

if res == 'no':
    particle_mass_max = float(input("Enter the maximum desired particle mass in solar masses:"))
    frame_duration = float(input("Enter the desired frame duration in milliseconds:"))
    t_final = float(input("Enter the time you would like the simulation to run for (myrs):"))
    L = float(input("Enter the length/height of the box (parsecs): "))/2
    num_bins = int(input("Enter the number of bins used in the correlation function: "))
    softening_length = float(input(" Enter the desired softening length (parsecs): "))
    Np = int(input("Enter the number of particles in the simulation: "))
    Nt = int(input("Enter the number of time steps used (sets deltat): "))
    Nd = int(input("Enter the number of spatial dimensions (2 or 3): "))
    v_max = float(input("Enter the maximum velocity of any of the particles (pcs/myrs): "))
    dist_type = input("Enter the type of distribution you would like the particles to initially have, "  
    "options are: 'u', 'g', 't' or 'f': " )
    n = 1
    if dist_type != 'u':
        n = float(input("Enter the standard deviation of any gaussian distribution (sets how big the " 
        "blobs of particles are): "))
    run_particle_sim(particle_mass_max = particle_mass_max, frame_duration = frame_duration, 
    t_final = t_final, L = L, num_bins = num_bins, softening_length = softening_length, Np = Np, 
    Nt = Nt, Nd = Nd, v_max = v_max, n = n, dist_type = dist_type)

if res == 'yes':
    run_particle_sim()



