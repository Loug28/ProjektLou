# Python molecular dynamics simulation of particles in 2 dimensions with real time animation
# BH, OF, MP 2020-10-20, latest verson 2021-10-24

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import random as rnd

"""
    This script is rather long: sit back and try to understand its structure before jumping into coding.
    MD simulations are performed by a class (MDsimulator) that envelops both the parameters and the algorithm;
    in this way, performing several MD simulations can be easily done by just allocating more MDsimulator
    objects instead of changing global variables and/or writing duplicates.
    You are asked to implement two things:
    - pair force and potential calculation and
    - temperature coupling.
    The latter is encapsulated into the class, so make sure you are modifying the variables and using the
    parameters of the class (the one you can access via 'self.variable_name' or 'self.function_name()').
"""

# Constant to obtain box dimensions given the number of particles per row 
crystal_spacing = 1.07

# Boltzmann constant
kB = 1.0

# Lennard-Jones constants
epsilon = 1
sigma_Lj = 1

# Number of steps between heat capacity output
N_OUTPUT_HEAT_CAP = 1000


# Generate two Gaussian random numbers with standard deviation sigma, mean 0
def gaussianRandomNumbers(sigma):
    w = 2
    while (w >= 1):
        rx1 = 2 * rnd.random() - 1
        rx2 = 2 * rnd.random() - 1
        w = rx1 * rx1 + rx2 * rx2
    w = math.sqrt(-2 * math.log(w) / w)
    return sigma * rx1 * w, sigma * rx2 * w


# Assigns Gaussian distributed velocities given an energy per particle
def thermalize(vx, vy, sqrtKineticEnergyPerParticle):
    for i in range(0, len(vx)):
        vx[i], vy[i] = gaussianRandomNumbers(sqrtKineticEnergyPerParticle)


# The pair potential
def pairEnergy(r):
    return 4 * epsilon * ((sigma_Lj / r) ** 12 - (sigma_Lj / r) ** 6)


# The pair force
def pairForce(r):
    inv_r = 1 / r
    return 24 * (2 * inv_r ** 13 - inv_r ** 7)


# Calculate the shortest periodic distance, unit cell [0,Lx],[0,Ly]
# Returns the difference along x, along y and the distance
# This code assumes all particles are within [0,Lx],[0,Ly]
def pbc_dist(x1, y1, x2, y2, Lx, Ly):
    dx = x1 - x2
    dy = y1 - y2
    while dx < -0.5 * Lx:
        dx += Lx
    while dx > 0.5 * Lx:
        dx -= Lx
    while dy < -0.5 * Ly:
        dy += Ly
    while dy > 0.5 * Ly:
        dy -= Ly
    return dx, dy, math.sqrt(dx * dx + dy * dy)


class MDsimulator:
    """
        This class encapsulates the whole MD simulation algorithm
    """

    def __init__(self,
                 n=24,
                 L=6.7,
                 mass=1.0,
                 numPerRow=6,
                 T=0.4,
                 dt=0.01,
                 nsteps=20000,
                 numStepsPerFrame=100,
                 startStepForAveraging=100
                 ):

        """
            This is the class 'constructor'; if you want to try different simulations with different parameters 
            (e.g. temperature) in the same scrip, allocate another simulator by passing a different value as input
            argument. See the examples at the end of the script.
        """

        # Initialize simulation parameters and box
        self.n = n
        self.mass = mass
        self.invmass = 1.0 / mass
        self.numPerRow = numPerRow
        self.Lx = L
        self.Ly = L
        self.T = T
        self.kBT = kB * T
        self.dt = dt
        self.nsteps = nsteps
        self.numStepsPerFrame = numStepsPerFrame
        # Initialize positions, velocities and forces
        self.x = []
        self.y = []
        self.vx = []
        self.vy = []
        self.fx = []
        self.fy = []
        for i in range(n):
            self.x.append(crystal_spacing * ((i % numPerRow) + 0.5 * (i / numPerRow)))
            self.y.append(crystal_spacing * 0.87 * (i / numPerRow))
            self.vx.append(0.0)
            self.vy.append(0.0)
            self.fx.append(0.0)
            self.fy.append(0.0)
        thermalize(self.vx, self.vy, np.sqrt(self.kBT / self.mass))
        # Initialize containers for energies
        self.sumEpot = 0
        self.sumEpot_squared = 0
        self.sumPV = 0
        self.outt = []
        self.ekinList = []
        self.epotList = []
        self.etotList = []
        self.pvList = []
        self.startStepForAveraging = startStepForAveraging
        self.step = 0
        self.Epot = 0
        self.Ekin = 0
        self.PV = 0
        self.Cv = 0
        self.use_andersen_thermostat = False
        # Initialize figure for animation
        # self.fig = plt.figure()
        # self.ax = plt.subplot(xlim=(0, self.Lx), ylim=(0, self.Ly))

    def clear_energy_potential(self):

        """
            Clear the temporary variables storing potential and kinetic energy
            Resets forces to zero
        """

        self.Epot = 0
        self.Ekin = 0
        for i in range(0, self.n):
            self.fx[i] = 0
            self.fy[i] = 0

    def update_forces(self):

        """
            Updates forces and potential energy using functions
            pairEnergy and pairForce (which you coded above...)
            Returns the virial: 1/2 sum f*r
        """

        sum_f_times_r = 0
        for i in range(self.n):
            for j in range(i + 1, self.n):
                dx, dy, r = pbc_dist(self.x[i], self.y[i], self.x[j], self.y[j], self.Lx, self.Ly)
                self.Epot += pairEnergy(r)
                fij = pairForce(r)
                self.fx[i] += fij * dx / r
                self.fy[i] += fij * dy / r
                self.fx[j] -= fij * dx / r
                self.fy[j] -= fij * dy / r
                sum_f_times_r += fij * r

        # Here we divide by 4 instead of 2 because we double-counted the forces
        return -sum_f_times_r / 4

    def propagate(self):

        """
            Performs an Hamiltonian propagation step and
            rescales velocities to match the input temperature 
            (THE LATTER YOU NEED TO IMPLEMENT!)
        """

        for i in range(0, self.n):
            # At the first step we already have the "full step" velocity
            if self.step > 0:
                # Update the velocities with a half step
                self.vx[i] += self.fx[i] * self.invmass * 0.5 * self.dt
                self.vy[i] += self.fy[i] * self.invmass * 0.5 * self.dt

            # When temperature coupling, modify the velocity of one particle here
            if self.use_andersen_thermostat and i % N_OUTPUT_HEAT_CAP == 0:
                thermalize(self.vx, self.vy, np.sqrt(self.kBT))

            # Add the kinetic energy of particle i to the total
            self.Ekin += 0.5 * self.mass * (self.vx[i] ** 2 + self.vy[i] ** 2)
            # Update the velocities with a half step
            self.vx[i] += self.fx[i] * self.invmass * 0.5 * self.dt
            self.vy[i] += self.fy[i] * self.invmass * 0.5 * self.dt
            # Update the coordinates
            self.x[i] += self.vx[i] * self.dt
            self.y[i] += self.vy[i] * self.dt
            # Apply p.c.b. and put particles back in the unit cell
            self.x[i] = self.x[i] % self.Lx
            self.y[i] = self.y[i] % self.Ly

    def md_step(self):

        """
            Performs a full MD step
            (computes forces, updates positions/velocities)
        """

        # This function performs one MD integration step
        self.clear_energy_potential()
        virial = self.update_forces()
        self.propagate()
        self.PV = 2 * (self.Ekin - virial)
        # Start averaging only after some initial spin-up time
        if self.step > self.startStepForAveraging:
            self.sumEpot += self.Epot
            self.sumEpot_squared += self.Epot * self.Epot
            self.sumPV += self.PV
        self.step += 1

    def integrate_some_steps(self, framenr=None):

        """
            Performs MD steps in a prescribed time window
            Stores energies and heat capacity
        """

        for j in range(self.numStepsPerFrame):
            self.md_step()
        t = self.step * self.dt
        self.outt.append(t)
        self.ekinList.append(self.Ekin)
        self.epotList.append(self.Epot)
        self.etotList.append(self.Epot + self.Ekin)
        self.pvList.append(self.PV)
        if self.step >= self.startStepForAveraging and self.step % N_OUTPUT_HEAT_CAP == 0:
            EpotAv = self.sumEpot / (self.step + 1 - self.startStepForAveraging)
            Epot2Av = self.sumEpot_squared / (self.step + 1 - self.startStepForAveraging)
            self.Cv = (Epot2Av - EpotAv * EpotAv) / (self.kBT * self.T)
            pvAv = self.sumPV / (self.step + 1 - self.startStepForAveraging)
            # print('time', t, '<Cv> =', self.Cv, "<P>V =", pvAv)

    def snapshot(self, framenr=None):

        """
            This is an 'auxillary' function needed by animation.FuncAnimation
            in order to show the animation of the 2D Lennard-Jones system
        """

        self.integrate_some_steps(framenr)
        return self.ax.scatter(self.x, self.y, s=1500, marker='o', c="r"),

    def simulate(self, use_andersen_thermostat=False):

        """
            Performs the whole MD simulation
            If the total number of steps is not divisible by the frame size, then
            the simulation will run for nsteps-(nsteps%numStepsPerFrame) steps
        """
        self.use_andersen_thermostat = use_andersen_thermostat

        nn = int(self.nsteps // self.numStepsPerFrame)
        # print("T=" + str(self.T) + ", integrating for " + str(nn * self.numStepsPerFrame) + " steps ...")
        for i in range(nn):
            self.integrate_some_steps()

    def simulate_animate(self):

        """
            Performs the whole MD simulation, while producing and showing the
            animation of the molecular system
            CAREFUL! This will slow down the script execution considerably
        """

        nn = self.nsteps // self.numStepsPerFrame
        print("T=" + str(self.T) + ", integrating for " + str(nn * self.numStepsPerFrame) + " steps ...")
        anim = animation.FuncAnimation(self.fig, self.snapshot,
                                       frames=nn, interval=50, blit=True, repeat=False)
        plt.show()  # show the animation
        # You may want to (un)comment the following 'waitforbuttonpress', depending on your environment
        # plt.waitforbuttonpress(timeout=20)

    def plot_energy(self, title='Simulation'):

        """
            Plots kinetic, potential and total energy over time
        """

        plt.figure()
        plt.xlabel('time')
        plt.ylabel('energy')
        plt.plot(self.outt, self.ekinList, self.outt, self.epotList, self.outt, self.etotList)  # , self.pvList)
        plt.legend(('Ekin', 'Epot', 'Etot'))  # , 'PV'))
        plt.title(title)
        plt.show()


def assignment_a():
    timesteps = [0.032, 0.034]
    simulations = {}
    for dt in timesteps:
        for i in range(5):
            sim = MDsimulator(dt=dt)
            sim.simulate()
            sim.plot_energy(title='Timestep dt = ' + str(dt))


def calculate_average(sims):
    """Input multiple list and returns a list with the mean"""
    avg_kinetic_energy = []
    avg_potential_energy = []
    avg_total_energy = []
    for sim in sims:
        avg_kinetic_energy.append(sim.ekinList)
        avg_potential_energy.append(sim.epotList)
        avg_total_energy.append(sim.etotList)
    result = {
        'avg_kinetic': np.average(np.array(avg_kinetic_energy), axis=0),
        'avg_potential': np.average(np.array(avg_potential_energy), axis=0),
        'avg_total': np.average(np.array(avg_total_energy), axis=0)
    }
    return result


def assignment_b():
    temperatures = [0.2, 1]
    simulations = {0.2: [], 1: []}
    result = {}
    for T in temperatures:
        for i in range(1):
            sim = MDsimulator(T=T)
            sim.simulate()
            simulations[T].append(sim)
        result[T] = calculate_average(simulations[T])
    plt.clf()
    plt.plot(simulations[0.2][0].outt, result[0.2]['avg_kinetic'], label='T=0.2')
    plt.plot(simulations[0.2][0].outt, result[1]['avg_kinetic'], label='T=1')
    plt.legend()
    plt.title('Kinetic energy over time for different temperatures')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.show()

    # <editor-fold desc="plots">
    plt.clf()
    plt.plot(simulations[0.2][0].outt, result[0.2]['avg_potential'], label='T=0.2')
    plt.plot(simulations[0.2][0].outt, result[1]['avg_potential'], label='T=1')
    plt.legend()
    plt.title('Potential energy over time for different temperatures')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.show()
    # </editor-fold>

    plt.clf()
    plt.plot(simulations[0.2][0].outt, result[0.2]['avg_total'], label='T=0.2')
    plt.plot(simulations[0.2][0].outt, result[1]['avg_total'], label='T=1')
    plt.legend()
    plt.title('Total energy over time for different temperatures')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.show()


def assignment_c():
    temperatures = [0.2, 1]
    simulations = {0.2: [], 1: []}
    result = {}
    for T in temperatures:
        for i in range(1):
            sim = MDsimulator(T=T)
            sim.simulate(use_andersen_thermostat=True)
            simulations[T].append(sim)
        result[T] = calculate_average(simulations[T])

    # <editor-fold desc="Plot">
    plt.clf()
    plt.plot(simulations[0.2][0].outt, result[0.2]['avg_kinetic'], label='T=0.2')
    plt.plot(simulations[0.2][0].outt, result[1]['avg_kinetic'], label='T=1')
    plt.legend()
    plt.title('Kinetic energy over time for different temperatures')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.show()

    plt.clf()
    plt.plot(simulations[0.2][0].outt, result[0.2]['avg_potential'], label='T=0.2')
    plt.plot(simulations[0.2][0].outt, result[1]['avg_potential'], label='T=1')
    plt.legend()
    plt.title('Potential energy over time for different temperatures')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.show()

    plt.clf()
    plt.plot(simulations[0.2][0].outt, result[0.2]['avg_total'], label='T=0.2')
    plt.plot(simulations[0.2][0].outt, result[1]['avg_total'], label='T=1')
    plt.legend()
    plt.title('Total energy over time for different temperatures')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.show()
    # </editor-fold>


def assignment_d():
    temperatures = [0.2, 0.3, 0.4, 0.45, 0.55, 0.6, 0.7, 0.75, 0.85, 1]
    simulations = {}
    result = {}
    avg_tot_energy = {}
    avg_heat_capacity = {}
    for T in temperatures:
        simulations[T] = []
        for i in range(2):
            sim = MDsimulator(T=T, startStepForAveraging=20000, nsteps=40000)
            sim.simulate(use_andersen_thermostat=True)
            simulations[T].append(sim)
        result[T] = calculate_average(simulations[T])
        avg_tot_energy[T] = np.mean(result[T]['avg_kinetic'])
        avg_heat_capacity[T] = np.mean([sim.Cv for sim in simulations[T]])

    plt.clf()
    plt.plot(temperatures, avg_tot_energy.values())
    plt.title('Average energy for simulation with temperature T')
    plt.xlabel('Temperatures')
    plt.ylabel('Average energy')
    plt.show()

    plt.clf()
    plt.plot(temperatures, avg_heat_capacity.values())
    plt.title('Heat capacity for simulation with temperature T')
    plt.xlabel('Temperatures')
    plt.ylabel('Heat capacity')
    plt.show()


def assignment_e():
    temperatures = np.array([0.2, 0.3, 0.4,0.5, 0.6,0.7, 0.8,0.9, 1])
    n = 24
    L = 6.7
    m = 2
    pressure_L_v_T = []
    for T in temperatures:
        pressures = np.zeros(m)
        for i in range(m):
            sim = MDsimulator(T=T)
            sim.simulate(use_andersen_thermostat=True)
            pressures[i] = sim.PV
        pressure_L_v_T.append(np.mean(pressures))  #pressure_L.append(sim.PV / ((2 * L) ** 3)) stod innan
    plt.clf()
    plt.plot(temperatures, pressure_L_v_T, label='Simulation')
    plt.plot(temperatures, n * kB * temperatures, label='Ideal gas')
    plt.legend()
    plt.xlabel('Temperatures')
    plt.ylabel('PV')
    plt.title('PV vs temperature for L = 6.7')
    plt.show()

    pressure_2L_v_T = []
    for T in temperatures:
        pressures = np.zeros(m)
        for i in range(m):
            sim = MDsimulator(T=T, L=2*L)
            sim.simulate(use_andersen_thermostat=True)
            pressures[i] = sim.PV
        pressure_2L_v_T.append(np.mean(pressures))
    plt.clf()
    plt.plot(temperatures, pressure_2L_v_T, label='Simulation')
    plt.plot(temperatures, n * kB * temperatures, label='Ideal gas')
    plt.legend()
    plt.xlabel('Temperatures')
    plt.ylabel('PV')
    plt.title('PV vs temperature for L = 13.4')
    plt.show()

# assignment_a()
# assignment_b()
# assignment_c()
# assignment_d()
assignment_e()