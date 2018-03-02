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

class Utility:
    def __init__(self):
        pass

    def get_fft_plot(self, a):
        """ INPUT: aperture object """
        screen_width = a.width
        fft_values = np.fft.fft(a.aperture_values)
        ks = np.linspace(0, screen_width, len(a.aperture_values)) # Use aperture dimension for scaling screen

        plt.plot(ks, np.fft.fftshift(fft_values)) # FFT reverses pattern through the middle
        plt.ylabel("Amplitude strength")
        plt.xlabel("metres")
        plt.show()

def main():
    a = Aperture(0.00001, 10)
    u = Utility()

    a.populate_aperturevalues_smallslit(0.0000001)
    u.get_fft_plot(a)

if __name__ == "__main__":
    main()