import numpy as np
import matplotlib.pyplot as plt
import random


class NVME:

    def __init__(self, L, N, T, dt, num_steps):
        # Define simulation parameters
        self.N = N  # number of particles
        self.Tbath = T  # temperature
        self.dt = dt  # time step
        self.rc = 2.5
        # Step 3: Generate initial conditions
        self.pos = np.zeros((N, 2))  # array to store particle positions
        self.lastpos = np.zeros((N, 2))  # array to store previous particles' positions
        self.vel = np.zeros((N, 2))  # array to store particle velocities
        self.forces = np.zeros((N, 2))  # array to store particle forces
        self.L = L  # size of the simulation domain

        # init for the radial distribution function g(r)
        self.CALCULATE_RDF = False
        self.r_array = []
        self.bin_width = 0.5  # width of each bin
        self.num_bins = int(self.L / 2 / self.bin_width)  # number of bins
        self.gr = np.zeros(self.num_bins)  # initialize g(r) to zeros
        self.counts = np.zeros(self.num_bins)  # initialize the counts for each bin to zeros

        self.Tins_list = []
        self.Tins = 0  # instantaneous temperature
        self.total_momentum = 0
        self.pe_list = []
        self.ke_list = []
        self.total_momentum_list = []
        self.num_steps = num_steps
        self.total_energy = 0
        self.ke = 0
        self.pe = 0
        self.total_energy_list = []  # list to store the total energy at each time step


        for i in range(N):
            self.pos[i] = [random.uniform(0, L),
                           random.uniform(0, L)]  # generate initial positions using uniform distribution
            self.vel[i] = np.random.normal(0, np.sqrt(self.Tbath),
                                           2)  # generate initial velocities using Maxwell-Boltzmann distribution
            # self.vel[i] = np.array([0, 0])

        self.vel -= np.mean(self.vel, axis=0)  # remove any net momentum

    # Define functions for Lennard-Jones potential and force
    def lj_potential(self, r, rc):
        if r <= rc:
            return (4 * ((1 / r) ** 12 - (1 / r) ** 6)) - (4 * ((1 / rc) ** 12 - (1 / rc) ** 6))
        else:
            return 0

    def lj_force(self, r, rc, direction):
        if r <= rc:
            f = ((48 * direction) / (r ** 2)) * ((1 / r) ** 12 - 0.5 * (1 / r) ** 6)
            return f
        else:
            return 0

    # Step 5: Implement the velocity-Verlet algorithm
    def verlet_integrate(self, pos, lastpos, heated_vel, forces):
        pe = 0  # for each time step make sure not to add previous steps potential energy
        ke = 0  # for each time step make sure not to add previous steps kinetic energy

        pos_new = 2 * pos - lastpos + 0.5 * forces * self.dt ** 2  # verlet position
        # vel_half = vel + 0.5 * forces * self.dt #frog algortihm
        vel_new = ((pos_new - lastpos) / 2 * self.dt)  # verlet velocity

        pos_new = np.mod(pos_new, self.L)  # apply periodic boundary conditions
        forces_new = np.zeros_like(forces)

        for i in range(self.N):
            for j in range(i + 1, self.N):
                r = pos_new[j] - pos_new[i]
                r = r - self.L * np.round(r / self.L)  # minimum image convention
                r_mag = np.linalg.norm(r)
                if r_mag != 0:
                    pe += self.lj_potential(r_mag, self.rc)  # potential
                    f = self.lj_force(r_mag, self.rc, r / r_mag)
                    forces_new[i] += f
                    forces_new[j] -= f
                    if self.CALCULATE_RDF:
                        # compute the bin index and increment the count for the corresponding bin
                        bin_index = int(r_mag / self.bin_width)
                        if bin_index < self.num_bins:
                            self.counts[bin_index] += 2 # counting particles
        # vel_new = vel_half + 0.5 * forces_new * self.dt #frog algorithm

        if self.CALCULATE_RDF:
            for i in range(self.num_bins):
                # compute the distance range for the bin
                r_min = i * self.bin_width
                r_max = (i + 1) * self.bin_width

                # compute the volume of the bin
                v_bin = (r_max ** 3 - r_min ** 3) * self.bin_width**3

                # compute the ideal gas g(r) value
                if i == 0 or self.counts[i] == 0:
                    gr_ideal = 1  # prevent division by zero
                else:
                    gr_ideal = (4 / 3) * np.pi * v_bin

                # compute the g(r) and r values
                self.gr[i] = (self.gr[i] + self.counts[i]) / (len(self.gr) * self.N * gr_ideal)
                self.r_array.append((r_max+r_min)/2)

        ke_current = np.sum(0.5 * np.sum(vel_new ** 2, axis=1))  # kinetic energy
        _Tins = ke_current * 2 / (3 * self.N)  # instantaneous temperature

        if self.Tbath > 0 and _Tins > 0:  # setting Tbath = 0 disables the thermostat
            # Berendsen thermostat
            # tau = self.dt / 0.0025
            scale_factor = 1 + (0.0025) * ((self.Tbath / _Tins) - 1)
            vel_new *= np.sqrt(scale_factor)
            ke_bath = np.sum(0.5 * np.sum(vel_new ** 2, axis=1))
            # Extract the magnitude of each velocity vector
            vel_mags = [np.linalg.norm(v) for v in vel_new]
            # Take the sum of the magnitudes to get the total momentum
            total_momentum = np.sum(vel_mags)

            total_energy = pe + ke_bath
            ke = ke_bath

        else:
            # Extract the magnitude of each velocity vector
            vel_mags = [np.linalg.norm(v) for v in vel_new]
            # Take the sum of the magnitudes to get the total momentum
            total_momentum = np.sum(vel_mags)
            # total_momentum = np.sum(vel_new, axis=0)

            total_energy = pe + ke_current
            ke = ke_current

        return pos_new, vel_new, forces_new, pe, ke, total_energy, total_momentum, _Tins, self.gr, self.r_array

    def simulate(self) -> None:
        # Run the simulation and plot the results
        for i in range(self.num_steps):
            print(i)
            if i == self.num_steps - 1:
                self.CALCULATE_RDF = True

            self.lastpos = self.pos
            self.pos, self.vel, self.forces, self.pe, self.ke, self.total_energy, self.total_momentum, self.Tins, self.gr, self.r_array = self.verlet_integrate(
                self.pos, self.lastpos, self.vel, self.forces)

            self.total_energy_list.append(self.total_energy)
            self.ke_list.append(self.ke)
            self.pe_list.append(self.pe)
            self.total_momentum_list.append(self.total_momentum)
            self.Tins_list.append(self.Tins)



    def plot_gr(self) -> None:
        plt.clf()
        plt.plot(self.r_array, self.gr)
        plt.xlabel('r')
        plt.ylabel('gr')
        plt.title(f"#steps={self.num_steps}, T={self.Tbath}, dt = {self.dt}, L = {self.L}, N = {self.N}")
        plt.savefig(f"gr#steps={self.num_steps}_T={self.Tbath}_dt = {self.dt}_L = {self.L}_N = {self.N}.png")
        # plt.show()

    def plot_Tins(self) -> None:
        plt.clf()
        plt.plot(np.arange(self.num_steps)[40:], self.Tins_list[40:])
        plt.xlabel('Time')
        plt.ylabel('Temperature')
        plt.title(f"#steps={self.num_steps}, T={self.Tbath}, dt = {self.dt}, L = {self.L}, N = {self.N}")
        plt.savefig(f"temperature#steps={self.num_steps}_T={self.Tbath}_dt = {self.dt}_L = {self.L}_N = {self.N}.png")
        # plt.show()

    def plot_total_energy(self) -> None:
        plt.clf()
        plt.plot(np.arange(self.num_steps)[40:], self.total_energy_list[40:])
        plt.xlabel('Time')
        plt.ylabel('Total energy')
        plt.title(f"#steps={self.num_steps}, T={self.Tbath}, dt = {self.dt}, L = {self.L}, N = {self.N}")
        plt.savefig(f"total_energy#steps={self.num_steps}_T={self.Tbath}_dt = {self.dt}_L = {self.L}_N = {self.N}.png")
        # plt.show()

    def plot_kinetic_energy(self) -> None:
        plt.clf()
        plt.plot(np.arange(self.num_steps)[40:], self.ke_list[40:])
        plt.xlabel('Time')
        plt.ylabel('KE')
        plt.title(f"#steps={self.num_steps}, T={self.Tbath}, dt = {self.dt}, L = {self.L}, N = {self.N}")
        plt.savefig(
            f"kinetic_energy#steps={self.num_steps}_T={self.Tbath}_dt = {self.dt}_L = {self.L}_N = {self.N}.png")
        # plt.show()

    def plot_potential_energy(self) -> None:
        plt.clf()
        plt.plot(np.arange(self.num_steps)[40:], self.pe_list[40:])
        plt.xlabel('Time')
        plt.ylabel('PE')
        plt.title(f"#steps={self.num_steps}, T={self.Tbath}, dt = {self.dt}, L = {self.L}, N = {self.N}")
        plt.savefig(
            f"potential_energy#steps={self.num_steps}_T={self.Tbath}_dt = {self.dt}_L = {self.L}_N = {self.N}.png")
        # plt.show()

    def plot_total_momentum(self) -> None:
        plt.clf()
        plt.plot(np.arange(self.num_steps)[40:], self.total_momentum_list[40:])
        plt.xlabel('Time')
        plt.ylabel('Total Momentum')
        plt.title(f"#steps={self.num_steps}, T={self.Tbath}, dt = {self.dt}, L = {self.L}, N = {self.N}")
        plt.savefig(
            f"total_momentum#steps={self.num_steps}_T={self.Tbath}_dt = {self.dt}_L = {self.L}_N = {self.N}.png")
        # plt.show()


# _numstepslist = [100, 625, 900, 2000, 100, 100, 625, 900, 2000]
# _T_list = [0, 0, 0, 0, 0.1, 1, 1, 1, 1]
# for i in range(len(_numstepslist)):
#     obj = NVME(30, 100, 1, 0.1, 500)
#     # obj = NVME(30, 100, _T_list[i], 0.1, _numstepslist[i])
#     obj.simulate()
#     # obj.plot_total_energy()
#     # obj.plot_potential_energy()
#     # obj.plot_kinetic_energy()
#     # obj.plot_total_momentum()
#     # obj.plot_Tins()
#     obj.plot_gr()
obj = NVME(30, 100, 1, 0.1, 500)
obj.simulate()
obj.plot_gr()

