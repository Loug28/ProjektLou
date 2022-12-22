#!/bin/python3

import random
import matplotlib.pyplot as plt
import numpy as np


class TrafficModel:
    def __init__(self, N, vmax, p, Ncars, roadLength, eqtime=0.2):
        self.currentN = 0
        self.N = N
        self.vmax = vmax
        self.p = p
        self.Ncars = Ncars
        self.roadLength = roadLength
        self.xList = list(range(0, Ncars))
        self.vList = Ncars * [0]
        self.flowRateList = []
        self.eqtime = eqtime
        self.run()

    def distance_to_next_car(self, car):
        # If we are the final car, we use the index of the first car.
        if car == len(self.xList) - 1:
            x_next = self.xList[0]
        else:
            x_next = self.xList[car + 1]

        return (x_next - self.xList[car] + self.roadLength) % self.roadLength

    def run(self):
        while (self.currentN != self.N):
            # self.plot() # plot every step
            # print(self.xList, self.vList) # debug

            for car in range(0, len(self.xList)):

                # Step 1
                if self.vList[car] < self.vmax:
                    self.vList[car] += 1

                # Step 2
                d = self.distance_to_next_car(car)
                if self.vList[car] >= d:
                    self.vList[car] = d - 1

                # Step 3
                if (random.random() < self.p) and (self.vList[car] >= 1):
                    self.vList[car] -= 1

                # Step 4
                self.xList[car] += self.vList[car]

            # Apply PBC
            for car in range(0, len(self.xList)):
                if self.xList[car] >= self.roadLength:
                    self.xList[car] -= self.roadLength

            # Record flowrate in list and discard first 20% for equilibration.
            if self.currentN > self.eqtime * self.N:
                self.flowRateList.append(sum(self.vList) / self.roadLength)

            # increment step
            self.currentN += 1

    # Allow the system to equilibrate before recording the flow rate.
    # Average the flow rate over multiple times to get more accurate results.
    def flowRate_mean(self):
        return sum(self.flowRateList) / len(self.flowRateList)

    def density(self):
        return self.Ncars / self.roadLength

    # Plot road on x-axis, time on y-axis and dots for cars.
    def plot(self):
        plt.scatter(self.Ncars * [self.currentN], self.xList)
        plt.ylim(0, 10)
        plt.savefig('test_{}.png'.format(self.currentN))
        plt.clf()
        plt.close()

################################################################################

def main_test():
    sim = TrafficModel(N=100, vmax=2, p=0.5, Ncars=5, roadLength=10)
    # sim.plot()

def main_2a():
    flowRateList = []
    densityList  = []

    for roadLength in range(5, 51):
        for Ncars in range(2, roadLength - 1):
            print(roadLength, Ncars)
            sim = TrafficModel(N=1000, vmax=2, p=0.5, Ncars=Ncars, roadLength=roadLength)

            flowRateList.append(sim.flowRate_mean())
            densityList.append(sim.density())
    
    plt.scatter(densityList, flowRateList)
    plt.xlabel('Density')
    plt.ylabel('Flowrate')
    plt.title('Fundametal diagram')
    plt.grid()
    plt.show()

def main_2b():
    flowRateList = []
    stdErrList = []
    OPointOOne = []
    
    while True:
        sim = TrafficModel(N=100, vmax=2, p=0.5, Ncars=25, roadLength=50, eqtime=0.2)
        flowRateList.append(sim.flowRate_mean())
        
        stderr = np.std(flowRateList) / len(flowRateList)**0.5
        stdErrList.append(stderr)
        OPointOOne.append(0.001)
        if stderr < 0.001 and len(flowRateList) > 2:
            print(len(flowRateList), stderr)
            break
    
    #plotting --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    plt.title("Statistical accuracy")
    plt.xlabel("Number of simulations")
    plt.ylabel("Standard errors of mean")
    plt.plot(range(0, len(stdErrList)), stdErrList, label="std error")
    plt.plot(range(0, len(stdErrList)), OPointOOne, label="0.001")
    plt.legend()
    plt.show()



def main_2c():
    flowRateList = []

    for roadLength in range(10, 1000):
        print
        sim = TrafficModel(N=1000, vmax=2, p=0.5, Ncars=20, roadLength=roadLength, eqtime=0.5)
        flowRateList.append(sim.flowRate_mean())
    
    plt.plot(range(10, 1000), flowRateList)
    plt.title('vmax = {}, p = {}, Ncars = {}, eqtime = {}'.format(sim.N, sim.vmax, sim.p, sim.eqtime))
    plt.xlabel('Road length')
    plt.ylabel('Flow rate')
    plt.show()

def main_2d():
    flowRateList1 = []
    flowRateList2 = []
    flowRateList3 = []
    densityList   = []

    for Ncars in range(2, 500):
        print(Ncars)
        
        sim1 = TrafficModel(N=100, vmax=1, p=0.5, Ncars=Ncars, roadLength=500, eqtime=0.3)
        sim2 = TrafficModel(N=100, vmax=2, p=0.5, Ncars=Ncars, roadLength=500, eqtime=0.3)
        sim3 = TrafficModel(N=200, vmax=5, p=0.5, Ncars=Ncars, roadLength=500, eqtime=0.3)

        flowRateList1.append(sim1.flowRate_mean())
        flowRateList2.append(sim2.flowRate_mean())
        flowRateList3.append(sim3.flowRate_mean())
        densityList.append(sim1.density())

    plt.plot(densityList, flowRateList1, label='vmax=1')
    plt.plot(densityList, flowRateList2, label='vmax=2')
    plt.plot(densityList, flowRateList3, label='vmax=5')
    plt.xlabel('Density')
    plt.ylabel('Flow rate')
    plt.legend()
    plt.show()

def main_2e():
    flowRateList1 = []
    flowRateList2 = []
    densityList   = []

    for roadLength in range(5, 51):
        for Ncars in range(2, roadLength - 1):
            print(roadLength, Ncars)
            sim1 = TrafficModel(N=5000, vmax=2, p=0.2, Ncars=Ncars, roadLength=roadLength, eqtime=0.5)
            sim2 = TrafficModel(N=5000, vmax=2, p=0.8, Ncars=Ncars, roadLength=roadLength, eqtime=0.5)

            flowRateList1.append(sim1.flowRate_mean())
            flowRateList2.append(sim2.flowRate_mean())
            densityList.append(sim1.density())
    
    plt.scatter(densityList, flowRateList1, label='p=0.2')
    plt.scatter(densityList, flowRateList2, label='p=0.8')
    plt.xlabel('Density')
    plt.ylabel('Flow rate')
    plt.title('Fundametal diagram')
    plt.legend()
    plt.grid()
    plt.show()

################################################################################

# main_test()
# main_2a()
#main_2b() # Need around 150 attempts for eqtime = 20%
main_2c()
# main_2d()
# main_2e()
