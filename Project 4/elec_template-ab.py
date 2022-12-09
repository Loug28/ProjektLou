import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.ndimage import gaussian_filter1d

L = 10

# The grid is n+1 points along x and y, including boundary points 0 and n
n = 20

# The grid spacing is L/n

# The number of iterations
nsteps = 10
nsteps_list = range(1,500)
error_max_array = []        # should have the same size as nsteps_list
v_exact = 10

# Initialize the grid to 0
# 10% lower potential 
v = np.ones((n+1, n+1)) * 9
vold = np.zeros((n+1, n+1)) * 9
vnew = np.zeros((n+1, n+1)) * 9

# 4.1b
#v = np.ones((n+1, n+1)) * 0

'''# 10.10b
for i in v:
    i[5] += 4'''

# Set the boundary conditions
for i in range(1,n):
    v[0,i] = 10
    v[n,i] = 10
    v[i,0] = 5
    v[i,n] = 5

#fig = plt.figure()
#ax = fig.add_subplot(111)
#im = ax.imshow(v, cmap=None, interpolation='nearest')
#fig.colorbar(im)

# checker=1: no checkboard, checker=2: checkerboard (note: n should be even)
checker = 2

# perform one step of relaxation
def relax(n, v, checker):
    
    error_array = []
    for check in range(0,checker):
        for x in range(1,n):
            for y in range(1,n):
                if (x*(n+1) + y) % checker == check:
                    v[x,y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1])*0.25
                    
                    error = abs(v_exact - v[x,y])
                    error_array.append(error)
        #error_max_array.append(max(error_array))

        # Copy back the new values to v
        # Note that you can directly store in v instead of vnew with Gauss-Seidel or checkerboard
        #for x in range(1,n):
        #    for y in range(1,n):
        #        if (x*(n+1) + y) % checker == check:
        #            v[x,y] = vnew[x,y]

        '''error_array = []
        for x in range(1,n):
            for y in range(1,n):
                if (x*(n+1) + y) % checker == check:
                    v[x,y] = vnew[x,y]

                    # Check 1% accuracy
                    error = abs(v[x,y]-vold[x,y])

                    average = (1/4) * (v[x+1, y] + v[x-1, y] + v[x, y+1] + v[x, y-1])
                    error = abs(average - v[x,y])
                    print("error:", error)
                    if error > maximumerror:
                        maximumerror = error
                        print("max:", maximumerror)
        error_max_array.append(maximumerror)'''
        

def update(step):
    #print(step)
    global n, v, checker

    # FuncAnimation calls update several times with step=0,
    # so we needs to skip the update with step=0 to get
    # the correct number of steps 
    if step > 0:
        relax(n, v, checker)

    #im.set_array(v)
    #return im,

def E4_1a():
    for step in nsteps_list:
        update(step)
    
    new = error_max_array[:len(error_max_array)//2]
    print(len(new))
    newsmooth = gaussian_filter1d(new, sigma=0.9)

    plt.figure(1)
    plt.title("Number of iterations to achive 1% accuracy")
    plt.xlabel("Number of iterations")
    plt.ylabel("Accuracy")
    plt.plot(nsteps_list, newsmooth, "g", label="Grid size: "+ str(L/n))
    plt.plot(nsteps_list, [0.1]*len(nsteps_list), 'r', label="1%")
    plt.legend()
    plt.show()


# we generate nsteps+1 frames, because frame=0 is skipped (see above)
#anim = animation.FuncAnimation(fig, update, frames=nsteps+1, interval=200, blit=True, repeat=False)
#plt.show()
E4_1a()

