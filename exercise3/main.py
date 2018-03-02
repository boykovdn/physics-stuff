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

    def populate_aperturevalues_smallslit(self, slitwidth):
        middle_index = len(self.aperture_values)//2
        M = np.int(slitwidth/self.delta)
        for i in range(0, M):
            self.aperture_values[i - int(M/2)] = 1

    def populate_aperturevalues_grating(self, gratingwidth):
        m = 8
        s = 100e-6

        fun = lambda x: np.exp((m / 2) * np.sin(2 * np.pi * x / s) * 1j)
        f_vect = np.vectorize(fun)

        vals = np.linspace(0, len(self.aperture_values), gratingwidth/self.delta)
        vals = f_vect(vals)
        print(vals.shape)

        middle_index = int(len(self.aperture_values)/2)
        self.aperture_values[middle_index - int(len(vals)/2) : middle_index - int(len(vals)/2) + len(vals)] = vals

class Utility:
    def __init__(self):
        pass

    def get_fft_plot(self, a):
        """ INPUT: aperture object """
        fft_values = np.fft.fft(a.aperture_values)

        # Make scaled x axis
        wavelength = 500e-9
        d = 100e-6
        D = 1
        L = 5e-3
        screen_width = a.width

        ns = np.linspace(0, len(a.aperture_values), len(a.aperture_values)) # Use aperture dimension for scaling screen
        xs_freq = np.fft.fftfreq(a.N, a.delta) * (wavelength * D)

        plt.plot(np.fft.fftshift(xs_freq), np.abs(np.fft.fftshift(fft_values))**2) # FFT reverses pattern through the middle
        plt.ylabel("Relative intensity (not normalised)")
        plt.xlabel("metres")
        plt.show()



def main():
    a = Aperture(1e-2, 10)
    u = Utility()

    a.populate_aperturevalues_grating(2e-3)
    u.get_fft_plot(a)

if __name__ == "__main__":
    main()