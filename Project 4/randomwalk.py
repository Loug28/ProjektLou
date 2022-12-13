import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit
from pylab import cm as plcm            # colors for plotting
import pickle                           # for storing calculations and reading again while editing plot details
import math
import random
from scipy.stats import sem
plt.rcParams.update({'font.size': 12})

# ------------------------- Set global values

n = 10
tolerance = 0.000001

Nlist = [10, 100, 250, 500, 1000, 2500, 5000, 7500, 10000]
x = 5
y = 6
V = 0

# ------------------------- Define functions

def initiateVMatrices():
    """Initiates potential matrixes with first boundary value: V = 10 on
    two opposing sides and 5 on the other two sides.
    - v has boundary values and initial guess of 9 everywhere else
    - vNew is a copy of v
    - vExact is the exact analytical solution, 10 everywhere"""
    global v, vNew
    # Initialize the grid to 0
    v = np.zeros((n+1, n+1))        # matrix of v, index are i: row, j:column
    # Set the boundary conditions
    for i in range(1,n):
        v[0,i] = 10
        v[n,i] = 10
        v[i,0] = 5
        v[i,n] = 5
    # Initial guess
    for i in range(1,n):
        for j in range(1,n):
            v[i,j] = 9
    vNew = np.copy(v)
    return vNew

def relax():
    """One relax iteration. v[i,j] is set as the avarage of its neighbours."""
    global v, vNew, n
    for x in range(1,n):
        for y in range(1,n):
            vNew[x,y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1])*0.25
    for x in range(1,n):
        for y in range(1,n):
            v[x,y] = vNew[x,y]

def randomwalk(x, y):
    global v, V
    for i in range(100):
        val = random.randint(0, 4) 
        if val == 0:
            x += 1

        if val == 1:
            y += 1

        if val == 2:
            x -= 1

        if val == 3:
            y -= 1
            
        if v[x,y] == 5 or v[x,y] == 10:
            V += v[x,y]
            break
        else:
            V += 0 
        

def calculate1():
    """Main calculation function that first initalizes with initiateVMatrixes()
    and then uses relax() until v is within tolerance.
    1. First boundary coditions
    2. Second boundary conditions"""
    global v, vNew, n, v1

    # First bondary conditions
    step = 0
    toleranceAcqurired = False
    while not toleranceAcqurired:
        step+=1
        vOld = np.copy(v)
        relax()
        # Controll accuracy
        toleranceAcqurired = True   # run through v and set false if not acquired
        for i in range(1,n):
            for j in range(1,n):
                if np.abs( (v[i,j]-vOld[i,j])/vOld[i,j] ) > tolerance:
                    toleranceAcqurired = False
    print('Tolerance was met after', step, 'steps.')
    v1 = np.copy(v)

# ----------------------- Plot

def calculate2():
    def walk(N):
        global V, v
        for i in range(N):
            randomwalk(x,y)    
        V1 = V
        V = 0
        return V1/N

    Vlist = []
    semlist = []
    errorlist = []
    exactlist = [v[x,y]]*len(Nlist)
    convergencelist = [walk(Nlist[-1])]*len(Nlist)
    for N in Nlist:
        Vlist.append(walk(N))
        sample = []
        for _ in range(10):       
            sample.append(walk(N))
        errorlist.append(sem(sample))
        

    plt.figure()
    plt.title('Convergence for (x,y) = (' 
              + str(x) + ', ' + str(y) + ')',
              fontsize = 18)
    plt.plot(Nlist, Vlist, color = 'r', markersize = 10)
    plt.plot(Nlist, exactlist, color = 'b', markersize = 10)
    plt.plot(Nlist, convergencelist, color = 'g', markersize = 10)
    plt.legend(['Random-walk solution',
                'Relaxation method convergence: V ≈ ' + str(round(v[x,y],1)),
                'Random walk convergence:  V ≈ ' + str(round(walk(Nlist[-1]),1))])
    plt.xlabel('Walkers (n)')
    plt.ylabel('Potential (V)')
    plt.show()

    plt.figure()
    plt.title('Standard error of the mean for (x,y) = (' 
              + str(x) + ', ' + str(y) + ')',
              fontsize = 18)
    plt.plot(Nlist, errorlist, marker = '.',
             color = 'tab:green', markersize = 10)
    plt.xlabel('Walkers (n)')
    plt.ylabel('SEM')
    plt.show()

initiateVMatrices()
calculate1()
calculate2()
print(v[x,y])