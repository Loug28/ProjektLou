import matplotlib.pyplot as plt
import numpy as np

class MonteCarlo:

    '''Class for initilizing Metropolis'''

    def __init__(self, delta, Nsteps, N0):
        self.delta = delta
        self.Nsteps = Nsteps
        self.N0 = N0
        self.x = []
        self.d = []