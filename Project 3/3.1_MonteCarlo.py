import matplotlib.pyplot as plt
import numpy as np
import math
import random as ran

class MonteCarlo:

    '''Class for initilizing Metropolis'''

    def __init__(self, delta, Nsteps, N0, dist):
        self.delta = delta
        self.Nsteps = Nsteps
        self.N0 = N0
        self.x = []
        self.d = []
        ran.randrange(0, 1)
        self.r = ran.random

    def P(x):
        if x > 0:
            return math.exp(-x)
        else:
            return 0

    def calc_x(self):
        xj = self.x.__getitem__(self.x.__len__ - 1) + self.d.__getitem__(self.d.__len__ - 1)
        if self.P(xj)/self.P(self.x.__getitem__(self.x.__len__ - 1)):
            self.x.append(xj)


""" import random

import matplotlib.pyplot as plt
import numpy as np


def e_to_the_power_of_negative_x(x):
    return np.exp(-x)


class MetropolisMethod:
    def __init__(self, delta, nsteps_total, nsteps_start_averaging, distribution):
        self.delta = delta
        self.nsteps_total = nsteps_total
        self.nsteps_start_averaging = nsteps_start_averaging
        self.distribution = distribution
        self.x_values = np.zeros(nsteps_total)

        self.calculate_x_values()

        self.result = np.sum(self.x_values[self.nsteps_start_averaging:])/(self.nsteps_total-self.nsteps_start_averaging)


    def calculate_x_values(self):
        for i in range(1, self.nsteps_total):
            d = random.uniform(-self.delta, self.delta)
            r = random.uniform(0, 1)
            if self.x_values[i-1] + d < 0 or self.distribution(self.x_values[i-1] + d) / self.distribution(self.x_values[i-1]) < r:
                self.x_values[i] = self.x_values[i - 1]
            else:
                self.x_values[i] = self.x_values[i-1] + d


def plot_result_v_delta():
    deltas = np.linspace(0.01, 1, 100)
    results = []
    for delta in deltas:
        results.append(MetropolisMethod(delta=delta, nsteps_total=30000, nsteps_start_averaging=1000,
                                        distribution=e_to_the_power_of_negative_x).result)
    plt.plot(deltas, results, '.-')
    plt.xlabel('delta')
    plt.ylabel('Calculated value on <x>')
    plt.title('<x> vs delta')
    plt.show()

def stat_error_plot():
    deltas = np.linspace(0.01, 10, 300)
    stat_error = {}
    for delta in deltas:
        stat_error[delta] = []
        for i in range(100):
            x = MetropolisMethod(delta=3, nsteps_total=20000, nsteps_start_averaging=1000,
                                distribution=e_to_the_power_of_negative_x).result
            stat_error[delta].append(x)
        stat_error[delta] = np.std(stat_error[delta])/10

    plt.clf()
    plt.plot(deltas, stat_error.values())
    plt.xlabel('delta')
    plt.ylabel('Statistical error')
    plt.title('Statistical error vs delta for 100 calculations')
    plt.show()

stat_error_plot()

 """