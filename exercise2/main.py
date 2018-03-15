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

        np.savetxt(filename,results_array,delimiter=",", header="time,angle,angular-velocity,energy")
        
class Plots:
    
    def __init__(self):
        self.u = Utility()

    def get_harmonic_angle(self, time):
        theta_0 = 0.01
        frequency = 1  # g=1, l=1
        return theta_0 * np.cos(frequency * time)

    def get_small_angle_oscillations(self):
        # Set G=0 (forcing). Second argument is damping, also 0
        # Get graphs for 10, 100, 1000 seconds
        ts10, coords10, energy10 = self.u.pendulum_damped([0.01, 0.0], 0, simulation_time=10*7, G=0)
        ts100, coords100, energy100 = self.u.pendulum_damped([0.01, 0.0], 0, simulation_time=100*7, G=0)
        ts1000, coords1000, energy1000 = self.u.pendulum_damped([0.01, 0.0], 0, simulation_time=1000*7, G=0)

        # Theoretical result for small-angle oscillations (no damping):
        # theta(t) = theta(0)*cos(wt)
        # theta(0) = theta_0 = 0.1

        time_theory = np.linspace(0, 1000*7, len(ts1000))
        theta_getter = np.vectorize(self.get_harmonic_angle)
        thetas = theta_getter(time_theory)

        f, (ax10, ax100, ax1000, ax_theory) = plt.subplots(4, sharex=True, sharey=True)
        ax10.plot(ts10, coords10[:, 0], label="10 periods sim")
        ax100.plot(ts100, coords100[:, 0], label="100 periods sim")
        ax1000.plot(ts1000, coords1000[:, 0], label="1000 periods sim")
        ax_theory.plot(time_theory, thetas, color="orange", label="Theoretical")
        ax10.legend(loc="best")
        ax100.legend(loc="best")
        ax1000.legend(loc="center right")
        ax_theory.legend(loc="center right")
        f.subplots_adjust(hspace=0)
        plt.show()
        
    def get_energy_plot(self):
        ts,coords,energy = self.u.pendulum_damped([0.01, 0.0], 0, simulation_time=100, G=0)
        plt.title("Energy evolution in time")
        plt.xlabel("time / s")
        plt.ylabel("energy / J")
        plt.plot(ts[:500],energy[:500])
        plt.text(ts[200], energy[50], "Error due to integrator.\nEnergy should be conserved otherwise")
        plt.show()
        # Error in rk4 algorithm

    # Period for Theta=pi/2 is roughly 3.7 seconds.
    def get_period_vs_amplitude(self):
        amplitudes = np.linspace(0.01, np.pi-0.001, 400) # Start from 0.01 because for 0 there is no motion
        velocities = []
        for a in amplitudes:
            ts, coords, energy = self.u.pendulum_damped([a,0.0], 0, simulation_time=40, G=0)
            velocities.append(coords[:, 1])

        period_indices = []
        print(period_indices)
        for v in velocities:
            period_indices.append(self.get_period_index(v))

        periods = []
        for i in period_indices:
            periods.append(ts[i])

        plt.scatter(amplitudes, periods, s=4)
        plt.ylabel("Period")
        plt.xlabel("Amplitude")
        plt.title("Period change as amplitude gets large")
        plt.text(x=0.1, y=14, s="Simulation stops at pi-0.001. Getting closer to\npi requires higher simulation times, since the pendulum\nis much slower to go down.")
        plt.show()

        #fft_angles = np.fft.fft(coords[:, 0])
        #fft_angles_abs = np.abs(fft_angles)[:int(len(fft_angles)/2)]
        #ts = ts[:int(len(ts)/2)]
        #plt.plot(ts, fft_angles_abs)
        #plt.show()

    def get_period_index(self, velocities):
        for i in range(1, len(velocities)-1):
            if velocities[i] * velocities[i+1] < 0:
                return i
        else:
            print("Problem with finding period from amplitude values")

    def damping_plot(self, dampings_list):
        for q in dampings_list:
            ts, coords, energy = self.u.pendulum_damped([0.5, 0], q, simulation_time=50, G=0)
            plt.plot(ts, coords[:, 1], label=("q="+str(q)))

        plt.xlabel("time")
        plt.ylabel("Angular velocity")
        plt.title("Damped oscillator")
        plt.legend()
        plt.show()

    def damping_forced_velocity(self, forcings_list):
        for f in forcings_list:
            ts, coords, energy = self.u.pendulum_damped([0.5,0], 0.5, G=f) # q = 0.5, constant
            plt.plot(ts, coords[:, 1], label=("F="+str(f)))

        plt.xlabel("time")
        plt.ylabel("Angular velocity")
        plt.title("Damped oscillator, different forcing")
        plt.text(x=0.2, y=-2.5, s="Period is constant. Response is a superposition\nof the transient and forced solution with"
                             "different periods.")
        plt.legend()
        plt.show()

    def damping_forced_angle(self, forcings_list):
        for f in forcings_list:
            ts, coords, energy = self.u.pendulum_damped([0.5, 0], 0.5, G=f)  # q = 0.5, constant
            plt.plot(ts, coords[:, 0], label=("F=" + str(f)))

        plt.xlabel("time")
        plt.ylabel("Angle")
        plt.title("Damped oscillator, different forcing")
        plt.legend()
        plt.show()

    def chaos_plot(self):
        t_sim = 10000
        ts0,coords0, energy0 = self.u.pendulum_damped([0.2, 0], 0, G=1.2, simulation_time=t_sim)
        ts1,coords1, energy1 = self.u.pendulum_damped([0.20001, 0], 0, G=1.2, simulation_time=t_sim,)
        plt.plot(ts0, coords0[:, 1], label="theta=0.2")
        plt.plot(ts1, coords1[:, 1], label="theta=0.20001")
        plt.xlabel("time / seconds")
        plt.ylabel("Angular velocity")
        plt.title("Initial conditions sensitivity test")
        plt.text(x=47.6, y=3.7, s="divergence point")
        plt.legend()
        plt.show()

    def chaos_angle_velocity(self):
        t_sim = 100
        ts, coords, energy = self.u.pendulum_damped([0, 0], 0, G=0.7, simulation_time=t_sim)
        #ts1, coords1, energy = self.u.pendulum_damped([0.001, 0], 0, G=0.001, simulation_time=t_sim)
        plt.plot(coords[:, 0], coords[:, 1])
        #plt.plot(coords1[:, 0], coords1[:, 1])
        plt.title("Phase space evolution of system")
        plt.ylabel("angle")
        plt.xlabel("angular velocity")
        plt.show()

def main():
    u = Utility()
    p = Plots()
    p.chaos_angle_velocity()

if __name__ == "__main__":
    main()

