import numpy as np
import matplotlib.pyplot as plt

class Aperture:
    """ Contains all functionality related to the aperture """
    def __init__(self, width, n):
        """ N=2**n where N is the number of samples of the aperture """
        self.N = 2**n
        self.width = width
        self.delta = self.width / self.N # The closest delta to delta_input that gives N as a power of 2
        self.aperture_values = np.zeros(self.N, dtype="complex")

    def set_aperture_values(self, complex_array):
        """ complex_array : 1d np array that holds the aperture values """
        if complex_array.shape != self.aperture_values.shape:
            print("Error: input array shape does not match aperture sample size (1,N)")
        else:
            self.aperture_values = complex_array

    def populate_aperturevalues_slit(self, slitwidth):
        """ Set the values in the self.aperture_values array to match a
        central slit that is slitwidth wide, with constant real amplitude. """
        M = np.int(slitwidth/self.delta)
        for i in range(-(M // 2), M // 2):
            self.aperture_values[self.N // 2 + i] = 1

    def populate_aperturevalues_grating(self, gratingwidth):
        m = 8
        s = 100e-6

        fun = lambda x: np.exp((m / 2) * np.sin(2 * np.pi * x / s) * 1j)
        f_vect = np.vectorize(fun)

        vals = np.linspace(0, len(self.aperture_values), gratingwidth/self.delta)
        vals = f_vect(vals)

        middle_index = int(len(self.aperture_values)/2)
        self.aperture_values[middle_index - int(len(vals)/2) : middle_index - int(len(vals)/2) + len(vals)] = vals

    def modify_nearfield(self, slitwidth, wavelength, D):
        """ Assume slit is centred in middle and modify values by a phase factor """

        k = (2 * np.pi) / wavelength

        fun = lambda x: np.exp(1j * (k * np.power(x, 2)) / (2 * D))
        fun_vect = np.vectorize(fun)

        coords = np.linspace(- self.width / 2, self.width / 2, self.N)
        phase_values = fun_vect(coords)
        self.aperture_values = self.aperture_values * phase_values

    def get_fft_plot(self, wavelength = 500e-9, D = 1):
        """ Return plot to be shown of diffraction pattern. Input D is the
        distance from the aperture to the screen. """
        fft_values = np.fft.fft(self.aperture_values)

        # Scale x axis
        xs = np.fft.fftfreq(self.N, self.delta) * (wavelength * D)

        plt.plot(np.fft.fftshift(xs), np.abs(np.fft.fftshift(fft_values))**2) # FFT reverses pattern through the middle
        plt.ylabel("Relative intensity (not normalised)")
        plt.xlabel("metres")

class Utility:
    """ Optional and helper functions. """

    def __init__(self):
        pass

    def get_fft_picture(self, a, filename):
        """ INPUT: aperture object """
        fft_values = np.fft.fft(a.aperture_values)

        # Make scaled x axis
        wavelength = 500e-9
        d = 100e-6
        D = 1
        L = 5e-3
        screen_width = a.width

        ns = np.linspace(0, len(a.aperture_values), len(a.aperture_values))  # Use aperture dimension for scaling screen
        xs_freq = np.fft.fftfreq(a.N, a.delta) * (wavelength * D)

        plt.plot(np.fft.fftshift(xs_freq),
                 np.abs(np.fft.fftshift(fft_values)) ** 2)  # FFT reverses pattern through the middle
        plt.ylabel("Relative intensity (not normalised)")
        plt.xlabel("metres")


    def compile_gif(self, dir):
        step_number = 0
        step_max = 80
        import imageio
        images = []
        while(step_number < step_max):
            images.append("../../practicals/gif_test/trippy_gif/" + "step" + str(step_number))
        imageio.mimsave(dir)

    def plot_theoretical(self, a, slitwidth):

        wavelength = 500e-9
        D = 1
        k = 2*np.pi/wavelength

        aperture_amplitude = 1
        scaling = 2*aperture_amplitude*(slitwidth/2)

        fun_theory = lambda x : scaling * np.sinc(x*(slitwidth/2)*(k/D))  # INPUT: distance across screen
        fun_theory_vect = np.vectorize(fun_theory)

        ys_n = np.linspace(-(a.N//2), a.N//2, a.N)
        y_length = ys_n * a.delta
        intensity = fun_theory_vect(y_length)**2

        visual_scaling = 1e23
        plt.plot(ys_n, visual_scaling*np.abs(intensity))

def phase_rbc(ys, l=5e-6, wavelength=500e-9):
    """
    l is the dimension of the rbc
    """
    n_rbc = 1.4 # Refractive index of red blood cell

    arg = lambda y : ((2*np.pi*n_rbc*l)/wavelength)*np.sin((np.pi/l)*y)**2
    exp_arg = lambda y : np.exp(1j * arg(y))
    vexp_arg = np.vectorize(exp_arg)
    phase = vexp_arg(ys)

    return phase

def phase_sphere(ys, l=5e-6, wavelength=500e-9):
    """
    l is the radius of the sphere
    """
    n_sphere = 1.4 # Refractive index of the sphere

    arg = lambda y : ((4*np.pi*n_sphere*l)/wavelength)*np.sqrt(1-(y/l)**2)
    exp_arg = lambda y : np.exp(1j * arg(y))
    vexp_arg = np.vectorize(exp_arg)
    phase = vexp_arg(ys)

    return phase

   


def main():
    #a = Aperture(1e-6, 20)
    #u = Utility()

    #slitwidth = 1e-7
    #a.populate_aperturevalues_slit(slitwidth)
    #u.plot_theoretical(a, slitwidth)
    #u.get_fft_plot(a)
    #plt.show()

    #a.populate_aperturevalues_slit(slitwidth)
    #a.modify_nearfield(slitwidth)
    #u.get_fft_plot(a)

    #counter = 0.00001
    #counter_name = 0
    #while counter < 0.0006:
    #    filename = "step" + str(counter_name)
    #    a.populate_aperturevalues_slit(counter)
    #    a.modify_nearfield(counter)
    #    u.get_fft_picture(a, "../../practicals/gif_test/" + filename)
    #    counter_name += 1
    #    counter += 0.000001

    #u.compile_gif("../../practicals/gif_test/trippy_gif/trippy_diffraction.gif")

    wavelength = 500e-9
    D = 2e-6 # Distance to imaging plane of microscope
    N = 19
    aperture_width = 1e-4

    """Core Task 1"""
    # L = 5mm
    # D = 1m
    # d = 100um
    # wavelength = 500nm

    # core1_pattern.png
    a = Aperture(aperture_width, N)

    l = 5e-6 # Dimension of rbc
    ys = np.linspace(-l,l,2**N)
    phase = phase_sphere(ys, l=l, wavelength=wavelength)

    a.set_aperture_values(phase) 
    #a.modify_nearfield(aperture_width, wavelength, D)
    a.get_fft_plot()
    plt.show()
    


if __name__ == "__main__":
    main()
