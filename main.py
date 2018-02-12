import scipy
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

class Utility():
    """Hold all functionality related to my integrator"""
    def __init__(self):
        self.g = 1

    def pendulum_undamped(self, coordinates, timestep=0.02):
        c = (9.81/1)
        func = lambda coords, ts_unused : [coords[1],-c*np.sin(coords[0])]
        #ts_unused because odeint requires the function to accept a second arg
        ts = np.linspace(0,10,10/timestep)
        ans = scipy.integrate.odeint(func,coordinates,ts)

        return ts,ans
    
    def pendulum_damped(self, coordinates, damping, timestep=0.02, length=1, mass=1, G=1, omega_D=1,simulation_time=50):
        q = damping/(mass*length)
        c = (self.g/length)
        F = G/(mass*np.power(length,2))
        func = lambda coords, t : [coords[1],(-c)*np.sin(coords[0]) + ((-q)*coords[1]) + (F*np.sin(omega_D*t))] #Returns [ddt_y0(t), ddt_y1(t)] to be passed to integrator
        ts = np.linspace(0,simulation_time,simulation_time/timestep)
        ans = scipy.integrate.odeint(func,coordinates,ts)
        print(ans)        

        energy = (np.power(ans[:,1],2)*0.5*mass*np.power(length,2)) + (length*mass*self.g*(1-np.cos(ans[:,0])))

        return ts,ans,energy
        
#Core Task 1: Writing out values
    def save_answers(self, filename, times, values,energy):
        """Export answers to "filename" """
        results_array = np.zeros((len(times),4))
        results_array[:,0] = times
        results_array[:,1] = values[:,0]
        results_array[:,2] = values[:,1]
        results_array[:,3] = energy

        np.savetxt(filename,results_array,delimiter=",")
        
class Plots:
    
    def __init__(self):
        self.u = Utility()

    def get_small_angle_oscillations(self, simulation_time):
        #Set G=0 (forcing). Second argument is damping, also 0
        ts,coords,energy = self.u.pendulum_damped([0.01,0.0],0,simulation_time=simulation_time,G=0)
        plt.plot(ts,coords[:,0])
        plt.show()
        

def main():
    p = Plots()
    p.get_small_angle_oscillations(3000)

if __name__ == "__main__":
    main()

