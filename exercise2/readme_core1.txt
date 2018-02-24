The relevant functions for this task are:

get_small_angle_oscillations()
get_energy_plot()
get_period_vs_amplitude()

which were used to generate the graphs.

All of these functions depend on:

pendulum_damped(self, coordinates, damping, timestep=0.02, length=1, mass=1, G=1, omega_D=1,simulation_time=50)

coordinates: [theta_initial, velocity_initial]

The other parameters are self-explanatory. I have separated the "back-end" computational functions in the Utility class
while the functions that generate plots have a separate Plots class. I think that helps code readability.

For theta_initial = pi/2 I get a period of roughly 3.7s