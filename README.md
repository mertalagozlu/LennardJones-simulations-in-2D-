# LennardJones-simulations-in-2D-
LennardJones simulations in 2D 


-Assume the simulation domain is a 2D square, L = 30 units in size, periodic in all directions.
- Use non-dimensional (reduced) units in your simulations (σ = 1, e = 1, rc=2.5) with truncated Lennard-Jones potential given as: 
The truncated Lennard-Jones potential is defined as:
Uₜᵣᵤₙcr= Ur- Urc,r ≤ rc; 0,r > rc
where U\left(r\right)=4\left[\left(\frac{1}{r}\right)^{12}-\left(\frac{1}{r}\right)^6\right], and r_c is the cutoff radius.
The directional-component of the Lennard-Jones force in the $x$-direction is given as:

F_x\left(r\right)=\frac{48x}{r^2}\left[\left(\frac{1}{r}\right)^{12}-0.5\left(\frac{1}{r}\right)^6\right]
and in the y-direction as:
F_y\left(r\right)=\frac{48y}{r^2}\left[\left(\frac{1}{r}\right)^{12}-0.5\left(\frac{1}{r}\right)^6\right]
- Assume that all particles have the same mass m = 1.
- Use velocity-Verlet scheme to integrate equations of motion in time.
- Use uniform distribution of N particles on the square lattice as your initial conditions. You can prescribe the initial particle velocities to be random, with the total momentum of the system equal to zero. Think about reasonable distribution of the velocities to be used if your initial temperature T should be close to some value (something between 0.1 and 1.0 may be a good value to start with). Note that the temperature can be estimated from the kinetic energy of the particles (note, that since we are using non-dimensional units in simulations, Boltzmann constant k_B should be taken equal to 1).
T=\frac{2}{3N}E_{kin}

where T is the temperature, N is the number of particles, and E_{kin} is the kinetic energy of the particles.
- Implementation and usage of cell list to speed up your simulations
- In molecular dynamics simulations, the total energy of a system can be calculated as the sum of its potential energy (PE) and kinetic energy (KE). The potential energy is the energy associated with the interactions between the particles in the system, while the kinetic energy is the energy associated with the motion of the particles.
The total energy formula for Lennard-Jones simulations in 3D can be as:
E_{\mathrm{total}}\ =\ E_{\mathrm{kinetic}}\ +\ E_{\mathrm{potential}}\ =\ \frac{1}{2}\ \sum_{i=1}^{N}{m_i\ {\mathbf{v}_\mathbf{i}}^\mathbf{2}}\ \ +\ 4\epsilon\ \sum_{i=1}^{N}\sum_{j>i}^{N}\left[\ \left(\frac{\sigma}{r_{ij}}\right)^{12}\ -\ \left(\frac{\sigma}{r_{ij}}\right)^6\ \right]\ \ 
where N is the number of particles, m_i and \mathbf{v}_\mathbf{i} are the mass and velocity of the ith particle, \epsilon and \sigma are parameters that determine the strength and distance of the interaction potential between particles, and r_{ij} is the distance between particles i and j. The first term on the right-hand side is the kinetic energy, and the second term is the potential energy.
- Implement calculation of radial distribution function. You can use the bin size of 0.05. If you decide to compute the radial distribution function for distances greater than cut-off radius rc, make sure that your link list can accommodate this (in case you combine computation of the radial distribution function with computation of pair-wise interactions).

The RDF is a metric for the likelihood of discovering a particle in a system at a given r-distance from another particle. The code must increase the count for each bin according to the separation of the particles in order to calculate the RDF.

The distance r_mag is divided by the self.bin_width to determine the bin index, which is then converted to an integer. Each bin's width in the RDF is specified by the self.bin_width attribute.
If the computed bin_index is less than self.num_bins, the code increases the count for the associated bin in the self.counts array by 2. The number 2 is chosen because the RDF is derived by counting pairs of particles. By increasing the count for the appropriate bin by 2, the code is essentially counting two particles separated by r.
- Implement Berendsen thermostat with relaxation time τ set as \Delta t/\tau = 0.0025.
scalefactor=\sqrt{\left(1+\Delta t/\tau\left(T_{\mathrm{bath}}/T_{\mathrm{ins}}-1\right)\right)}

The class NVME is defined, which takes in simulation parameters such as the size of the simulation domain, the number of particles, the temperature of the bath, the time step size, and the number of simulation steps to run.

By creating initial conditions for the particle locations and velocities, the __init__ function configures the simulation. 

To guarantee that the system is at rest, the velocities are taken from a Maxwell-Boltzmann distribution after the net momentum has been subtracted. The potential and force functions of Lennard-Jones are described. The velocity-Verlet algorithm is implemented by the verlet_integrate function, which modifies the particle locations and velocities at each time step. The Lennard-Jones potential and force functions are used to calculate the forces, while the Verlet method is used to update the locations and velocities. The positions are subjected to periodic boundary conditions. At each time step, the system's kinetic energy, potential energy, instantaneous temperature, and total energy are all determined. The velocities are rescaled in accordance with the Berendsen thermostat if a thermostat is enabled. Additionally, the system's overall momentum is calculated. Particle locations, velocities, kinetic and potential energies, temperature in the moment, total energy, and total momentum are all output by the simulation at each time step.

