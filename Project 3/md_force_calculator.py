# Python molecular dynamics simulation of particles in 2 dimensions with real time animation
# BH, OF, MP, AJ, TS 2022-11-20, latest verson 2021-10-21

import numpy as np
import random as rnd
import math

# Using numba to speed up force calculation
# More info: https://numba.pydata.org/numba-doc/latest/user/5minguide.html
import numba as nmb

"""
    Implement potential and force calculation
"""
LJ_epsilon = 1
LJ_sigma = 1

@nmb.jit(nopython=True)
def gaussianRandomNumbers(sigma):
    w = 2
    while (w >= 1):
        rx1 = 2 * rnd.random() - 1
        rx2 = 2 * rnd.random() - 1
        w = rx1 * rx1 + rx2 * rx2 
    w = np.sqrt(-2 * np.log(w) / w)
    return sigma * rx1 * w, sigma * rx2 * w

# Assigns random velocity components to all particles taken from a Gaussian
# distribution with sigma = sqrtKineticEnergyPerParticle, mu = 0
@nmb.jit(nopython=True)
def thermalize(vx, vy, sqrtKineticEnergyPerParticle):
    for i in range(0, len(vx)):
        vx[i], vy[i] = gaussianRandomNumbers(sqrtKineticEnergyPerParticle)

# The pair potential
@nmb.jit(nopython=True)
def pairEnergy(r):
    # Implement the pair energy calculation here!
    V = 4*LJ_epsilon*((LJ_sigma/r)**12-(LJ_sigma/r)**6)
    return V

# The pair force
@nmb.jit(nopython=True)    
def pairForce(r):
    # Implement the pair force calculation here!
    #F = -4*LJ_epsilon*((-12/r)*((LJ_sigma/r)**12) + (6/r)*((LJ_sigma/r)**6))
    inv_r = 1 / r
    F = 24 * (2 * inv_r ** 13 - inv_r ** 7)
    return F

# Calculate the shortest periodic distance, unit cell [0,Lx],[0,Ly]
# Returns the difference along x, along y and the distance
# This code assumes all particles are within [0,Lx],[0,Ly]
@nmb.jit(nopython=True)
def pbc_dist(x1, y1, x2, y2, Lx, Ly):
    dx = x1 - x2
    dy = y1 - y2
    while dx < -0.5*Lx:
        dx += Lx
    while dx > 0.5*Lx:
        dx -= Lx
    while dy < -0.5*Ly:
        dy += Ly
    while dy > 0.5*Ly:
        dy -= Ly
    return dx, dy, math.sqrt(dx*dx + dy*dy)

@nmb.jit(nopython=True)
def quick_force_calculation(x, y, fx, fy, Lx, Ly, n) :
    Epot = 0.0
    Virial = 0.0
    for i in range(n):
            for j in range(i+1,n):
                dx, dy, r = pbc_dist(x[i], y[i], x[j], y[j], Lx, Ly)
                Epot += pairEnergy(r)
                fij = pairForce(r)
                Virial -= 0.5*fij*r
                fx[i] += fij * dx / r
                fy[i] += fij * dy / r
                fx[j] -= fij * dx / r
                fy[j] -= fij * dy / r
    return Epot, Virial