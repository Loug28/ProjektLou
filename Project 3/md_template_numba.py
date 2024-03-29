# Python molecular dynamics simulation of particles in 2 dimensions with real time animation
# BH, OF, MP, AJ, TS 2022-11-20, latest verson 2021-10-21

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import random as rnd

# This local library contains the functions needed to perform force calculation
# Since this is by far the most expensive part of the code, it is 'wrapped aside'
# and accelerated using numba (https://numba.pydata.org/numba-doc/latest/user/5minguide.html)
import md_force_calculator as md

"""

    This script is rather long: sit back and try to understand its structure before jumping into coding.
    MD simulations are performed by a class (MDsimulator) that envelops both the parameters and the algorithm;
    in this way, performing several MD simulations can be easily done by just allocating more MDsimulator
    objects instead of changing global variables and/or writing duplicates.

    You are asked to implement two things:
    - Pair force and potential calculation (in md_force_calculator.py)
    - Temperature coupling (in md_template_numba.py)
    The latter is encapsulated into the class, so make sure you are modifying the variables and using the
    parameters of the class (the one you can access via 'self.variable_name' or 'self.function_name()').

"""

# Boltzmann constant
kB = 1.0

# Number of steps between heat capacity output
N_OUTPUT_HEAT_CAP = 1000

# You can use this global variable to define the number of steps between two applications of the thermostat
N_STEPS_THERMO = 10

# Lower (increase) this if the size of the disc is too large (small) when running run_animate()
DISK_SIZE = 750

class MDsimulator:

    """
        This class encapsulates the whole MD simulation algorithm
    """

    def __init__(self, 
        n = 48, 
        mass = 1.0, 
        numPerRow = 8, 
        initial_spacing = 1.12,
        T = 1, 
        dt = 0.01, 
        nsteps = 20000, 
        numStepsPerFrame = 100,
        startStepForAveraging = 100
        ):
        
        """
            This is the class 'constructor'; if you want to try different simulations with different parameters 
            (e.g. temperature, initial particle spacing) in the same scrip, allocate another simulator by passing 
            a different value as input argument. See the examples at the end of the script.
        """

        # Initialize simulation parameters and box
        self.n = n
        self.mass = 1.0
        self.invmass = 1.0/mass
        self.numPerRow = numPerRow
        self.Lx = numPerRow*initial_spacing
        self.Ly = numPerRow*initial_spacing
        self.area = self.Lx*self.Ly
        self.T = T
        self.kBT = kB*T
        self.dt = dt
        self.nsteps = nsteps
        self.numStepsPerFrame = numStepsPerFrame
        # Initialize positions, velocities and forces
        self.x = []
        self.y = []
        for i in range (n):
            self.x.append(self.Lx*0.95/numPerRow*((i % numPerRow) + 0.5*(i/numPerRow)))
            self.y.append(self.Lx*0.95/numPerRow*0.87*(i/numPerRow))
        
        # Numba likes numpy arrays much more than list
        # Numpy arrays are mutable, so can be passed 'by reference' to quick_force_calculation
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.vx = np.zeros(n, dtype=float)
        self.vy = np.zeros(n, dtype=float)
        self.fx = np.zeros(n, dtype=float)
        self.fy = np.zeros(n, dtype=float)

        # Initialize particles' velocity according to the initial temperature
        md.thermalize(self.vx, self.vy, np.sqrt(self.kBT/self.mass))
        # Initialize containers for energies
        self.sumEkin = 0
        self.sumEpot = 0
        self.sumEtot = 0
        self.sumEtot2 = 0
        self.sumVirial = 0
        self.outt = []
        self.ekinList = []
        self.epotList = []
        self.etotList = []
        self.startStepForAveraging = startStepForAveraging
        self.step = 0
        self.Epot = 0
        self.Ekin = 0
        self.Virial = 0
        self.Cv = 0
        self.P = 0
        self.use_andersen_thermostat = False

    def clear_energy_potential(self) :
        
        """
            Clear the temporary variables storing potential and kinetic energy
            Resets forces to zero
        """
        
        self.Epot = 0
        self.Virial = 0
        for i in range(0, self.n):
            self.fx[i] = 0
            self.fy[i] = 0

    def update_forces(self) :

        """
            Updates forces and potential energy using functions
            pairEnergy and pairForce (which you coded above...)
        """
        
        tEpot, tVirial = md.quick_force_calculation(self.x, self.y, self.fx, self.fy, 
            self.Lx, self.Ly, self.n)
        self.Epot += tEpot
        self.Virial += tVirial
    
    def propagate(self) :

        """
            Performs an Hamiltonian propagation step and
            rescales velocities to match the input temperature 
            (THE LATTER YOU NEED TO IMPLEMENT!)
        """

        # TODO
        # When using a thermostat, modify the velocities of all particles here.
        # Note that you can use thermalize() from md_force_calculator.py.
        #md.thermalize(self.vx, self.vy, np.sqrt(self.kBT))
        self.Ekin = 0
        for i in range(0,self.n):
            # At the first step we alread have the "full step" velocity
            if self.step > 0:
                # Update the velocities with a half step
                self.vx[i] += self.fx[i]*self.invmass*0.5*self.dt
                self.vy[i] += self.fy[i]*self.invmass*0.5*self.dt

            if self.use_andersen_thermostat and i % N_STEPS_THERMO == 0:
                md.thermalize(self.vx, self.vy, np.sqrt(self.kBT))
            # Add the kinetic energy of particle i to the total
            self.Ekin += 0.5*self.mass*(self.vx[i]*self.vx[i] + self.vy[i]*self.vy[i])
            # Update the velocities with a half step
            self.vx[i] += self.fx[i]*self.invmass*0.5*self.dt
            self.vy[i] += self.fy[i]*self.invmass*0.5*self.dt
            # Update the coordinates
            self.x[i] += self.vx[i] * self.dt
            self.y[i] += self.vy[i] * self.dt
            # Apply p.c.b. and put particles back in the unit cell
            self.x[i] = self.x[i] % self.Lx
            self.y[i] = self.y[i] % self.Ly

    def md_step(self) :

        """
            Performs a full MD step
            (computes forces, updates positions/velocities)
        """

        # This function performs one MD integration step
        self.clear_energy_potential()
        self.update_forces()
        # Start averaging only after some initial spin-up time
        if self.step > self.startStepForAveraging:
            self.sumVirial += self.Virial
            self.sumEkin   += self.Ekin
            self.sumEpot   += self.Epot
            self.sumEtot   += self.Epot+self.Ekin
            self.sumEtot2  += (self.Epot+self.Ekin)*(self.Epot+self.Ekin)
        self.propagate()
        self.step += 1

    def integrate_some_steps(self, framenr=None) :

        """
            Performs MD steps in a prescribed time window
            Stores energies and heat capacity
        """

        for j in range(self.numStepsPerFrame) :
            self.md_step()
        t = self.step*self.dt
        self.outt.append(t)
        self.ekinList.append(self.Ekin)
        self.epotList.append(self.Epot)
        self.etotList.append(self.Epot + self.Ekin)
        if self.step >= self.startStepForAveraging and self.step % N_OUTPUT_HEAT_CAP == 0:
            EkinAv  = self.sumEkin/(self.step + 1 - self.startStepForAveraging)
            EtotAv = self.sumEtot/(self.step + 1 - self.startStepForAveraging)
            Etot2Av = self.sumEtot2/(self.step + 1 - self.startStepForAveraging)
            VirialAV = self.sumVirial/(self.step + 1 - self.startStepForAveraging)
            self.Cv = (Etot2Av - EtotAv * EtotAv) / (self.kBT * self.T)
            self.P = (2.0/self.area)*(EkinAv - VirialAV)
            print('time', t, 'Cv =', self.Cv, 'P = ', self.P)

    def snapshot(self, framenr=None) :

        """
            This is an 'auxillary' function needed by animation.FuncAnimation
            in order to show the animation of the 2D Lennard-Jones system
        """

        self.integrate_some_steps(framenr)
        return self.ax.scatter(self.x, self.y, s=DISK_SIZE, marker='o', c="r"),

    def simulate(self, use_andersen_thermostat = False) :

        """
            Performs the whole MD simulation
            If the total number of steps is not divisible by the frame size, then
            the simulation will undergo nsteps-(nsteps%numStepsPerFrame) steps
        """
        self.use_andersen_thermostat = use_andersen_thermostat
        nn = self.nsteps//self.numStepsPerFrame
        print("Integrating for "+str(nn*self.numStepsPerFrame)+" steps...")
        for i in range(nn) :
            self.integrate_some_steps()

    def simulate_animate(self) :

        """
            Performs the whole MD simulation, while producing and showing the
            animation of the molecular system
            CAREFUL! This will slow down the script execution considerably
        """

        self.fig = plt.figure()
        self.ax = plt.subplot(xlim=(0, self.Lx), ylim=(0, self.Ly))

        nn = self.nsteps//self.numStepsPerFrame
        print("Integrating for "+str(nn*self.numStepsPerFrame)+" steps...") 
        self.anim = animation.FuncAnimation(self.fig, self.snapshot,
            frames=nn, interval=50, blit=True, repeat=False)
        plt.axis('square')
        plt.show()  # show the animation
        # You may want to (un)comment the following 'waitforbuttonpress', depending on your environment
        # plt.waitforbuttonpress(timeout=20)

    def plot_energy(self, title=("Energies")) :
        
        """
            Plots kinetic, potential and total energy over time
        """
        
        plt.figure(title)
        plt.title(title)
        plt.xlabel('time')
        plt.ylabel('energy')
        plt.plot(self.outt, self.ekinList, self.outt, self.epotList, self.outt, self.etotList)
        plt.legend( ('Ekin','Epot','Etot') )
        plt.savefig(title + ".pdf")
        plt.show()


# It's good practice to encapsulate the script execution in 
# a main() function (e.g. for profiling reasons)
def exercise_32a() :
    timesteps = np.arange(0.010, 0.031, 0.005)
    simulations = {}
    for dt in timesteps:
        for i in range(2):
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

def exercise_32b():
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

def exercise_32c():
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

def exercise_32d():
    temperatures = [0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8, 0.9, 1]
    simulations = {}
    result = {}
    avg_tot_energy = {}
    avg_heat_capacity = {}
    for T in temperatures:
        simulations[T] = []
        for i in range(2):
            sim = MDsimulator(T=T, startStepForAveraging=1000, nsteps=50000)
            sim.simulate(use_andersen_thermostat=True)
            simulations[T].append(sim)
        result[T] = calculate_average(simulations[T])
        avg_tot_energy[T] = np.mean(result[T]['avg_total'])
        avg_heat_capacity[T] = np.mean([sim.Cv for sim in simulations[T]])
        sim.simulate_animate()
    
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

def exercise_32e():
    temperatures = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    initial_spacing = 1.12
    m = 10 #10
    pressure_L_v_T = []
    for T in temperatures:
        pressures = np.zeros(m)
        for i in range(m):
            sim = MDsimulator(T=T)
            sim.simulate(use_andersen_thermostat=True)
            pressures[i] = sim.P
        pressure_L_v_T.append(np.mean(pressures))  #pressure_L.append(sim.PV / ((2 * L) ** 3)) stod innan
    plt.clf()
    plt.plot(temperatures, pressure_L_v_T, label='Simulation')
    plt.plot(temperatures, sim.n * kB * temperatures, label='Ideal gas')
    plt.legend()
    plt.xlabel('Temperatures')
    plt.ylabel('P')
    plt.title('P vs temperature')
    plt.show()

    pressure_2L_v_T = []
    for T in temperatures:
        pressures = np.zeros(m)
        for i in range(m):
            sim = MDsimulator(T=T, initial_spacing = (2 * initial_spacing))
            sim.simulate(use_andersen_thermostat=True)
            pressures[i] = sim.P
        pressure_2L_v_T.append(np.mean(pressures))
    plt.clf()
    plt.plot(temperatures, pressure_2L_v_T, label='Simulation')
    plt.plot(temperatures, sim.n * kB * temperatures, label='Ideal gas')
    plt.legend()
    plt.xlabel('Temperatures')
    plt.ylabel('P')
    plt.title('P vs temperature for doubling L')
    plt.show()

def main():
    #exercise_32a()
    #exercise_32b()
    #exercise_32c()
    exercise_32d()
    #exercise_32e()

# Calling 'main()' if the script is executed.
# If the script is instead just imported, main is not called (this can be useful if you want to
# write another script importing and utilizing the functions and classes defined in this one)
if __name__ == "__main__" :
    main()
