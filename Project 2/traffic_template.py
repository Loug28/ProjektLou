#!/bin/python3

# Template for traffic simulation
# BH, MP 2021-11-15, latest version 2022-11-1.

"""
    This template is used as backbone for the traffic simulations.
    Its structure resembles the one of the pendulum project, that is you have:
    (a) a class containing the state of the system and it's parameters
    (b) a class storing the observables that you want then to plot
    (c) a class that propagates the state in time (which in this case is discrete), and
    (d) a class that encapsulates the aforementioned ones and performs the actual simulation
    You are asked to implement the propagation rule(s) corresponding to the traffic model(s) of the project.
"""

import math
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy.random as rng
import numpy as np

import matplotlib


class Cars:

    """ Class for the state of a number of cars """

    def __init__(self, numCars=5, roadLength=50, v0=1):
        self.numCars    = numCars
        self.roadLength = roadLength
        self.t  = 0
        self.x  = []
        self.v  = []
        self.c  = []
        for i in range(numCars):
            self.x.append(i)        # the position of the cars on the road
            self.v.append(v0)       # the speed of the cars
            self.c.append(i)        # the color of the cars (for drawing)

    def distance(self, i):
        # Om bil x2 - x1 > 0 så måste man kolla på distansen framför mellan X2 och X1 D2 (Eftersom periodisk)
        # Om bil x2-x1 < 0 så måste man kolla på distansen mellan X1 och X2 D1
        # TODO: Implement the function returning the PERIODIC distance 
        # between car i and the one in front
        distance = (self.x[(i+1)%self.numCars] - self.x[i]) % self.roadLength
        return distance



class Observables:

    """ Class for storing observables """

    def __init__(self):
        self.time = []          # list to store time
        self.flowrate = []      # list to store the flow rate
        

class BasePropagator:

    def __init__(self):
        return
        
    def propagate(self, cars, obs):

        """ Perform a single integration step """
        
        fr = self.timestep(cars, obs)

        # Append observables to their lists
        obs.time.append(cars.t)
        obs.flowrate.append(fr)  # CHANGE!
              
    def timestep(self, cars, obs):

        """ Virtual method: implemented by the child classes """
        
        pass
      
        
class ConstantPropagator(BasePropagator) :
    
    """ 
        Cars do not interact: each position is just 
        updated using the corresponding velocity 
    """
    
    def timestep(self, cars, obs):
        for i in range(cars.numCars):
            cars.x[i] += cars.v[i]
        cars.t += 1
        return 0

# TODO
# HERE YOU SHOULD IMPLEMENT THE DIFFERENT CAR BEHAVIOR RULES
# Define you own class which inherits from BasePropagator (e.g. MyPropagator(BasePropagator))
# and implement timestep according to the rule described in the project

class MyPropagator(BasePropagator) :

    def __init__(self, vmax, p):
        BasePropagator.__init__(self)
        self.vmax = vmax
        self.p = p

    def timestep(self, cars, obs):
        # TODO Here you should implement the car behaviour rules
        for i in range(cars.numCars):
            if cars.v[i] < self.vmax:
                cars.v[i] = cars.v[i] + 1

        for j in range(cars.numCars):
            if cars.v[j] >= cars.distance(j):
                cars.v[j] = cars.distance(j) - 1

        for k in range(cars.numCars):
            if np.random.rand() <= self.p  and cars.v[k] > 0:
                cars.v[k] = cars.v[k] - 1

        for l in range(cars.numCars):
            cars.x[l] += cars.v[l]
            cars.x[l] %= cars.roadLength

        flowrate = np.sum(cars.v) / cars.roadLength
        cars.t += 1
        return flowrate



############################################################################################

def draw_cars(cars, cars_drawing):

    """ Used later on to generate the animation """
    theta = []
    r     = []

    for position in cars.x:
        # Convert to radians for plotting  only (do not use radians for the simulation!)
        theta.append(position * 2 * math.pi / cars.roadLength)
        #theta.append(position / cars.roadLength)
        r.append(1)

    return cars_drawing.scatter(theta, r, c=cars.c, cmap='hsv')


def animate(framenr, cars, obs, propagator, road_drawing, stepsperframe):

    """ Animation function which integrates a few steps and return a drawing """

    for it in range(stepsperframe):
        propagator.propagate(cars, obs)

    return draw_cars(cars, road_drawing),


class Simulation:

    def reset(self, cars=Cars()) :
        self.cars = cars
        self.obs = Observables()

    def __init__(self, cars=Cars()) :
        self.reset(cars)

    def plot_observables(self, title="simulation"):
        plt.clf()
        plt.title(title)
        plt.plot(self.obs.time, self.obs.flowrate)
        plt.xlabel('time')
        plt.ylabel('flow rate')
        plt.savefig(title + ".pdf")
        plt.show()

    # Run without displaying any animation (fast)
    def run(self,
            propagator,
            numsteps=200,           # final time
            title="simulation",     # Name of output file and title shown at the top
            ):

        for it in range(numsteps):
            propagator.propagate(self.cars, self.obs)

        self.plot_observables(title)

    # Run while displaying the animation of bunch of cars going in circe (slow-ish)
    def run_animate(self,
            propagator,
            numsteps=200,           # Final time 
            stepsperframe=1,        # How many integration steps between visualising frames
            title="simulation",     # Name of output file and title shown at the top
            ):

        numframes = int(numsteps / stepsperframe)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='polar')
        ax.axis('off')
        # Call the animator, blit=False means re-draw everything
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[self.cars,self.obs,propagator,ax,stepsperframe],
                                       frames=numframes, interval=50, blit=True, repeat=False)
        plt.show()

        # If you experience problems visualizing the animation and/or
        # the following figures comment out the next line 
        #plt.waitforbuttonpress(30)

        self.plot_observables(title)


#Assignment 2.2 a)
########################################################################################
flowratearr = []
densityarr = []
numOfSim = 1000
numOfSteps = 100
roadlen = 50
rangeStart = 0
rangeUpperLimit = 50
rangeStepSize = 2

def fund_calc():
    for numOfCars in range(rangeStart, rangeUpperLimit, rangeStepSize):
        flowrateArr_numOfCars = []
        for sim in range(numOfSim):
            cars = Cars(numCars=numOfCars, roadLength=roadlen)
            obs = Observables()
            propagator = MyPropagator(vmax=2, p=0.5)
            for it in range(numOfSteps):
                propagator.propagate(cars, obs)
            flowrateArr_numOfCars.append(np.sum(obs.flowrate)/numOfSteps)
        frcheck = np.sum(flowrateArr_numOfCars)/numOfSim
        print(str(frcheck) + "\n")
        flowratearr.append(frcheck)
        densityarr.append(numOfCars/roadlen)
        print ("Done with simulation for " + str (numOfCars) + " number of cars")

def plot_fund():
    plt.clf
    plt.title("Fundamental diagram")
    plt.xlabel("Density")
    plt.ylabel("Flow rate")
    plt.plot(densityarr, flowratearr) 
    plt.savefig("Fundamental.pdf")
    plt.show()


#Assignment 2.2 b)
######################################################################
meanFlowrate = []
error = []
numIter = 1

def stat_acc_round():
    global numIter
    while True:
        frArr = []
        for it in range(1, numIter, 1):
            cars = Cars(numCars=25, roadLength=50)
            obs = Observables()
            propagator = MyPropagator(vmax=2, p=0.5)
            for it in range(100):
                propagator.propagate(cars, obs)
            frArr.append(obs.flowrate)
        mean = np.mean(frArr)
        print (str(mean) + "\n")
        sem = mean / np.sqrt(np.size(frArr))
        print (str(sem) + "\n")
        error.append(sem)
        numIter += 1
        if round(sem, 3) <= 0.001:
            break
        
def stat_acc():
    global numIter
    while True:
        frArr = []
        for it in range(1, numIter, 1):
            cars = Cars(numCars=25, roadLength=50)
            obs = Observables()
            propagator = MyPropagator(vmax=2, p=0.5)
            for it in range(100):
                propagator.propagate(cars, obs)
            frArr.append(obs.flowrate)
        mean = np.mean(frArr)
        print (str(mean) + "\n")
        sem = mean / np.sqrt(np.size(frArr))
        print (str(sem) + "\n")
        error.append(sem)
        numIter += 1
        if sem <= 0.001:
            break

def stat_plot():
    plt.clf
    plt.title("Statistical accuracy")
    plt.xlabel("Number of simulations")
    plt.ylabel("Standard errors of mean")
    plt.plot(range(1, numIter, 1), error)
    plt.show()
    


#Assignment 2.2 c)
#######################################################################

roadLenArr = []
numSteps = 100
rStart = 10
rUpperLimit = 150
rIncrease = 20
densityUpper = 100
denArr = []

def independent_check():

    
    for length in range(rStart, rUpperLimit, rIncrease):
        denArr.clear()
        frLengthArr = []
        for den in range(0, densityUpper, 1):
            density = den / densityUpper
            savedfr = []
            for sim in range(numOfSim):
                cars = Cars(numCars=round(length*density), roadLength=length)
                obs = Observables()
                propagator = MyPropagator(vmax=2, p=0.5)
                for it in range(numSteps):
                    propagator.propagate(cars, obs)
                savedfr.append(np.sum(obs.flowrate)/numOfSteps)
            frcheck = np.sum(savedfr) / numOfSim
            print(str(frcheck) + "\n")
            frLengthArr.append(frcheck)
            denArr.append(density)
            print(str(density) + "\n")
        roadLenArr.append(frLengthArr)


            


def independent_plot():
    plt.clf()
    plt.xlabel("Density")
    plt.ylabel("Flow rate")
    plt.title("Fundamental diagram")
    for it in range(len(roadLenArr)):
        plt.plot(denArr, roadLenArr[it])
    plt.show()

# It's good practice to encapsulate the script execution in 
# a main() function (e.g. for profiling reasons)
def main() :

    # Here you can define one or more instances of cars, with possibly different parameters, 
    # and pass them to the simulator 

    # Be sure you are passing the correct initial conditions!
    cars = Cars(numCars = 10, roadLength=50)

    # Create the simulation object for your cars instance:
    simulation = Simulation(cars)

    #simulation.run_animate(propagator=ConstantPropagator())
    simulation.run_animate(propagator=MyPropagator(vmax=2, p=0.5))

#Assignment 2.2 a)
    #fund_calc()
    #plot_fund()

#Assignment 2.2 b)
    #stat_acc_round()
    #stat_plot()

#Assignment 2.2 c)
    independent_check()
    independent_plot()






# Calling 'main()' if the script is executed.
# If the script is instead just imported, main is not called (this can be useful if you want to
# write another script importing and utilizing the functions and classes defined in this one)
if __name__ == "__main__" :
    main()

