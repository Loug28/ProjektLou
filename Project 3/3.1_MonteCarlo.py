import matplotlib.pyplot as plt
import numpy as np
import math
import random as ran

class MonteCarlo:

    '''Class for initilizing Metropolis'''

    def __init__(self, delta, Nsteps, N0):
        self.delta = delta
        self.Nsteps = Nsteps
        self.N0 = N0
        self.x = np.zeros(Nsteps)
        self.d = []

        self.r = ran.uniform(0.0, 1.0)
        self.result = np.sum(self.x[N0:]) / (Nsteps - N0)
        """ print(str(self.result)) """

    def P(self, x):
        if x >= 0:
            return math.exp(-x)
        else:
            return 0

    def calc_x(self):
        for i in range(1, self.Nsteps):
            xi = self.x.__getitem__(i - 1)
            di = ran.uniform(-self.delta, self.delta)
            xj = xi + di
            if (self.P(x = xj) / self.P(x = xi)) > self.r:
                self.x[i] = xj
                """ print(str(xj)) """
            else:
                self.x[i] = (xi)
        self.result = np.sum(self.x[self.N0:]) / (self.Nsteps - self.N0)

steps = 200
def result_v_deltas():
    deltas = np.linspace(0.01, 10, steps)
    results = []
    for delta in deltas:
        x = MonteCarlo(delta=delta, Nsteps=10000, N0=1000)
        x.calc_x()
        results.append(x.result)
    plt.plot(deltas, results, '.-')
    plt.xlabel('delta')
    plt.ylabel('<x>')
    plt.title('<x> compared to delta')
    plt.show()

def SE_of_the_mean():
    deltas = np.linspace(0.01, 10, steps)
    errors = np.zeros(steps)
    ActualDiff = np.zeros(steps)
    io = 0
    for delta in deltas:
        print(str(delta))
        error = []
        for i in range(100):
            x = MonteCarlo(delta=delta, Nsteps=10000, N0=1000)
            x.calc_x()
            error.append(x.result)
        errors[io] = np.std(error) / 10 #10 is the square root of 100 and as we use N = 100 here with the range we will divide by 10
        ActualDiff[io] = errors[io] ** 0.5 / delta
        io += 1

    #plot the results
    #######################################################################'
    plt.clf()
    plt.plot(deltas, errors, label = "sigma/sqrt(N)")
    plt.plot(deltas, ActualDiff, label = "sqrt(sigma)/delta")
    plt.legend()
    plt.xlabel('delta')
    plt.ylabel('Exact value')
    #plt.title('SE vs delta')
    plt.show()


def main():
    #result_v_deltas()
    SE_of_the_mean()

if __name__ == "__main__" :
    main()