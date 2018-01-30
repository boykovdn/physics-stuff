import scipy
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np


class Utility():
    """Hold all functionality related to my integrator"""
    def __init__(self):
        pass

    def pendulum_undamped(self, coordinates, timestep=0.02):
        c = (9.81/1)
        func = lambda coords, ts_unused : [coords[1],-c*np.sin(coords[0])]
        #ts_unused because odeint requires the function to accept a second arg
        ts = np.linspace(0,10,10/timestep)
        ans = scipy.integrate.odeint(func,coordinates,ts)

        return ts,ans
    
    def pendulum_damped(self, coordinates, damping, timestep=0.02, length=1, mass=1, G=1, omega_D=1,simulation_time=50):
        q = damping/(mass*length)
        c = (9.81/length)
        F = G/(mass*np.power(length,2))
        func = lambda coords, t : [coords[1],-c*np.sin(coords[0]) - q*coords[1] + F*np.sin(omega_D*t)]
        ts = np.linspace(0,simulation_time,simulation_time/timestep)
        ans = scipy.integrate.odeint(func,coordinates,ts)

        return ts,ans
        

def main():
    u = Utility()
    ts,ans = u.pendulum_damped([0.02,0],timestep=0.002,damping=0.1,simulation_time=300,G=0) 
    plt.plot(ts,ans)
    plt.show()

if __name__ == "__main__":
    main()
