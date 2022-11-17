#!/bin/python3

# Python simulation of a simple planar pendulum with real time animation
# BH, OF, MP, AJ, TS 2020-10-20, latest version 2022-10-25.

from matplotlib import animation
from pylab import *
import copy

"""
    This script defines all the classes needed to simulate (and animate) a single pendulum.
    Hierarchy (somehow in order of encapsulation):
    - Oscillator: a struct that stores the parameters of an oscillator (harmonic or pendulum)
    - Observable: a struct that stores the oscillator's coordinates and energy values over time
    - BaseSystem: harmonic oscillators and pendolums are distinguished only by the expression of
                    the return force. This base class defines a virtual force method, which is
                    specified by its child classes
                    -> Harmonic: specifies the return force as -k*t (i.e. spring)
                    -> Pendulum: specifies the return force as -k*sin(t)
    - BaseIntegrator: parent class for all time-marching schemes; function integrate performs
                    a numerical integration steps and updates the quantity of the system provided
                    as input; function timestep wraps the numerical scheme itself and it's not
                    directly implemented by BaseIntegrator, you need to implement it in his child
                    classes (names are self-explanatory)
                    -> EulerCromerIntegrator: ...
                    -> VerletIntegrator: ...
                    -> RK4Integrator: ...
    - Simulation: this last class encapsulates the whole simulation procedure; functions are 
                    self-explanatory; you can decide whether to just run the simulation or to
                    run while also producing an animation: the latter option is slower
"""

# Global constants
G = 9.8  # gravitational acceleration

class Oscillator:

    """ Class for a general, simple oscillator """

    def __init__(self, m=1, c=4, t0=0, theta0=0.1 * np.pi, dtheta0=0, gamma=0):
        self.m = m              # mass of the pendulum bob
        self.c = c              # c = g/L
        self.L = G / c          # string length
        self.t = t0             # the time
        self.theta = theta0     # the position/angle
        self.dtheta = dtheta0   # the velocity
        self.gamma = gamma      # damping

class Observables:

    """ Class for storing observables for an oscillator """

    def __init__(self):
        self.time = []          # list to store time
        self.pos = []           # list to store positions
        self.vel = []           # list to store velocities
        self.energy = []        # list to store energy


class BaseSystem:
    
    def force(self, osc):

        """ Virtual method: implemented by the childc lasses  """

        pass


class Harmonic(BaseSystem):
    def force(self, osc):
        return - osc.m * ( osc.c*osc.theta + osc.gamma*osc.dtheta )


class Pendulum(BaseSystem):
    def force(self, osc):
        return - osc.m * ( osc.c*np.sin(osc.theta) + osc.gamma*osc.dtheta )


class BaseIntegrator:

    def __init__(self, _dt=0.01) :
        self.dt = _dt   # time step

    def integrate(self, simsystem, osc, obs):

        """ Perform a single integration step """
        
        self.timestep(simsystem, osc, obs)

        # Append observables to their lists
        obs.time.append(osc.t)
        obs.pos.append(osc.theta)
        obs.vel.append(osc.dtheta)
        # Function 'isinstance' is used to check if the instance of the system object is 'Harmonic' or 'Pendulum'
        if isinstance(simsystem, Harmonic) :
            # Harmonic oscillator energy
            obs.energy.append(0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + 0.5 * osc.m * G * osc.L * osc.theta ** 2)
        else :
            # Pendulum energy
            # TODO: Append the total energy for the pendulum (use the correct formula!)
            obs.energy.append(0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + osc.m * G * osc.L * (1 - cos(osc.theta)))


    def timestep(self, simsystem, osc, obs):

        """ Virtual method: implemented by the child classes """
        
        pass


# HERE YOU ARE ASKED TO IMPLEMENT THE NUMERICAL TIME-MARCHING SCHEMES:

class EulerCromerIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m
        osc.t += self.dt

        # TODO: Implement the integration here, updating osc.theta and osc.dtheta
        osc.theta = osc.theta + osc.dtheta * self.dt
        osc.dtheta = osc.dtheta - osc.c * osc.theta * self.dt


class VerletIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m
        osc.t += self.dt

        # TODO: Implement the integration here, updating osc.theta and osc.dtheta

        osc.theta = osc.theta + osc.dtheta * self.dt + 0.5 * accel * self.dt ** 2
        accel1 = simsystem.force(osc) / osc.m
        osc.dtheta = osc.dtheta + 0.5 * (accel1 + accel) * self.dt
    


class RK4Integrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m 
        # osc.t += self.dt

        # TODO: Implement the integration here, updating osc.theta and osc.dtheta
        
        theta_temp = copy.deepcopy(osc.theta)
        dtheta_temp = copy.deepcopy(osc.dtheta)
        t_temp = copy.deepcopy(osc.t)

        a1 = accel * self.dt
        b1 = dtheta_temp * self.dt

        osc.theta = theta_temp + (b1 * 0.5)
        osc.dtheta = dtheta_temp + (a1 * 0.5)
        osc.t = t_temp + (self.dt * 0.5)

        a2 = (simsystem.force(osc) / osc.m) * self.dt
        b2 = (dtheta_temp + (a1 * 0.5)) * self.dt

        osc.theta = theta_temp + (b2 * 0.5)
        osc.dtheta = dtheta_temp + (a2 * 0.5)

        a3 = (simsystem.force(osc)  / osc.m) * self.dt
        b3 = (dtheta_temp + (a2 * 0.5)) * self.dt

        osc.theta = theta_temp + b3
        osc.dtheta = dtheta_temp + a3
        osc.t = t_temp + self.dt

        a4 = (simsystem.force(osc) / osc.m) * self.dt
        b4 = (dtheta_temp + a3) * self.dt

        osc.dtheta = dtheta_temp + ((a1 + (2 * a2) + (2 * a3) + a4) / 6)
        osc.theta = theta_temp + ((b1 + (2 * b2) + (2 * b3) + b4) / 6)




# Animation function which integrates a few steps and return a line for the pendulum
def animate(framenr, simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe):
    
    for it in range(stepsperframe):
        integrator.integrate(simsystem, oscillator, obs)

    x = np.array([0, np.sin(oscillator.theta)])
    y = np.array([0, -np.cos(oscillator.theta)])
    pendulum_line.set_data(x, y)
    return pendulum_line,


class Simulation:

    def reset(self, osc=Oscillator()) :
        self.oscillator = osc
        self.obs = Observables()

    def __init__(self, osc=Oscillator()) :
        self.reset(osc)

    # Run without displaying any animation (fast)
    def run(self,
            simsystem,
            integrator,
            tmax=30.,               # final time
            ):

        n = int(tmax / integrator.dt)

        for it in range(n):
            integrator.integrate(simsystem, self.oscillator, self.obs)

    # Run while displaying the animation of a pendulum swinging back and forth (slow-ish)
    # If too slow, try to increase stepsperframe
    def run_animate(self,
            simsystem,
            integrator,
            tmax=30.,               # final time
            stepsperframe=5         # how many integration steps between visualising frames
            ):

        numframes = int(tmax / (stepsperframe * integrator.dt))-2

        # WARNING! If you experience problems visualizing the animation try to comment/uncomment this line
        plt.clf()

        # If you experience problems visualizing the animation try to comment/uncomment this line
        # fig = plt.figure()

        ax = plt.subplot(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_line, = ax.plot([], [], lw=5)
        plt.title(title)
        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[simsystem,self.oscillator,self.obs,integrator,pendulum_line,stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)

        # If you experience problems visualizing the animation try to comment/uncomment this line
        # plt.show()

        # If you experience problems visualizing the animation try to comment/uncomment this line
        plt.waitforbuttonpress(10)

    # Plot coordinates and energies (to be called after running)
    def plot_observables(self, title="simulation", ref_E=None):

        plt.clf()
        plt.title(title)
        plt.plot(self.obs.time, self.obs.pos, 'b-', label="Position")
        plt.plot(self.obs.time, self.obs.vel, 'r-', label="Velocity")
        plt.plot(self.obs.time, self.obs.energy, 'g-', label="Energy")
        if ref_E != None :
            plt.plot([self.obs.time[0],self.obs.time[-1]] , [ref_E, ref_E], 'k--', label="Ref.")
        plt.xlabel('time')
        plt.ylabel('observables')
        plt.legend()
        plt.savefig(title + ".pdf")
        plt.show()


# It's good practice to encapsulate the script execution in 
# a function (e.g. for profiling reasons)
def exercise_11() :
    # TODO
    sim11 = Simulation()

    sim11.run(simsystem=Pendulum(), integrator=VerletIntegrator())
    sim11.plot_observables(title = "Pendulum Verlet, gamma = 0")


def exercise_12():
    t_pendulum = []
    t_harmonic = []
    t_pertubation = []
    theta0_12 = np.linspace(0.1*np.pi, 0.9*np.pi, 1500)

    integrator = VerletIntegrator()

    for i in theta0_12:
        osc12 = Oscillator(theta0=i)
        obs12 = Observables()

        while np.sign(osc12.theta)==1:
            integrator.integrate(simsystem=Harmonic(), osc=osc12, obs=obs12)
        t_harmonic.append(osc12.t*4)


        osc12_P = Oscillator(theta0=i)

        while np.sign(osc12_P.theta)==1:
            integrator.integrate(simsystem=Pendulum(), osc=osc12_P, obs=obs12)
        t_pendulum.append(osc12_P.t*4)

        t_pertubation.append(2 * np.pi * np.sqrt(osc12.L / G) * (1 + (1/16 * i ** 2) + (11/3072 * i ** 4) + (173/737280 * i ** 6)))

    plt.figure()
    plt.plot(theta0_12, t_pendulum)
    plt.plot(theta0_12, t_harmonic)
    plt.plot(theta0_12, t_pertubation)
    plt.title("Perturbation series comparision to Pendulum and Harmonic")
    plt.figlegend(('Pendulum', 'Harmonic', 'Perturbation to theta(0) ^ 6'))
    plt.xlabel('Position')
    plt.ylabel('Time (s)')
    plt.show()


 
def exercise_13():
    #TODO
    sim13 = Simulation()

    sim13.run(simsystem=Harmonic(), integrator=VerletIntegrator())
    sim13.plot_observables(title = "Harmonic Verlet, gamma = 3")

def exercise_14():
    #TODO
    osc14 = Oscillator()
    sim14 = Simulation()
    obs14 = Observables()
    integrate = VerletIntegrator()
    
    time14 = []
    theta14 = []
    dtheta14 = []
    
    
    while osc14.t < 30:
        theta14.append(osc14.theta)
        dtheta14.append(osc14.dtheta)
        time14.append(osc14.t)
        integrate.integrate(simsystem=Pendulum(), osc=osc14, obs=obs14)
      
    plt.figure()
    plt.plot(dtheta14, theta14)

    plt.title("Phase space portrait of the position and velocity")
    plt.xlabel('Velocity')
    plt.ylabel('Position')
    plt.show()


    """ 
        This directive instructs Python to run what comes after ' if __name__ == "__main__" : '
        if the script pendulum_template.py is executed 
        (e.g. by running "python3 pendulum_template.py" in your favourite terminal).
        Otherwise, if pendulum_template.py is imported as a library 
        (e.g. by calling "import pendulum_template as dp" in another Python script),
        the following is ignored.
        In this way you can choose whether to code the solution to the exericises here in this script 
        or to have (a) separate script(s) that include pendulum_template.py as library.
    """

if __name__ == "__main__" :
    #exercise_11()
    exercise_12()
    #exercise_13()
    #exercise_14()
