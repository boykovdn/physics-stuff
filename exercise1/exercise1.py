import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


# N random points in d-dimensional space, values in [0,a)
def rng_points(a, d, N):
    result = np.random.rand(N, d)
    result = a * result
    return result


# 8-dimensional sine function sin(x0,x1,...,x7)
def sine_8dim(vector):
    return (np.sin(np.sum(vector)))


# Vectorized extraction of function values, apply function() on the rows of the matrix function_arg_array and return the vector Fi = function(Xi).
def get_function_values(function, function_arg_array):
    result = np.apply_along_axis(function, 1, function_arg_array)
    return result


def get_results_to_plot(function, s, d, N):
    V = np.power(s, d)  # d-dimensional cube volume
    f_vals = get_function_values(function, rng_points(s, d, N))

    avg = lambda array: np.sum(array) / array.shape[0]
    squarer = lambda t: t * t

    square = np.vectorize(squarer)
    f_vals_squared = square(f_vals)

    f_avg = avg(f_vals)
    f_avg_sq = np.power(f_avg, 2)
    f_sq_avg = avg(f_vals_squared)

    integral = f_avg * V * np.power(10, 6)
    sigma = np.sqrt((f_sq_avg - f_avg_sq) / f_vals.shape[0]) * V * np.power(10, 6)

    return integral, sigma


# higher data_number seems to give a better result
def get_sigma_stat(function, s, d, N, data_number=50):
    V = np.power(s, d)
    integral_values = []

    for i in range(0, data_number):
        f_vals = get_function_values(function, rng_points(s, d, N))
        f_avg = np.average(f_vals)
        integral_values.append(f_avg * V * np.power(10, 6))

    return np.std(integral_values)


# def graph_fit(sigmas,Ns):
#    slope,intercept,R_value,P_value,std_err = stats.linregress(sigmas,Ns)
#    return slope,intercept

# Supplementary Task 1
def get_std_comparison_plot():
    s = np.pi / 8
    fs = []
    sigmas_theoretical = []
    sigmas_stat = []
    N = 1000
    Ns = []

    for i in range(0, 25):
        Ns.append(N)
        f, sigma_theoretical = get_results_to_plot(sine_8dim, s, 8, N)
        sigma_stat = get_sigma_stat(sine_8dim, s, 8, N)
        print(sigma_theoretical, sigma_stat)
        fs.append(f)
        sigmas_theoretical.append(sigma_theoretical)
        sigmas_stat.append(sigma_stat)
        N += 500

    plt.plot(Ns, sigmas_theoretical, label="Theoretical stdev")
    plt.plot(Ns, sigmas_stat, label="Statistical stdev")
    plt.legend(loc=1)
    plt.xlabel("N, number of samples")
    plt.show()


# Core Task 2
def get_fresnel_integrals(u, data_points=100):
    us = np.linspace(-u, u, data_points)

    fresnel_from_zero = lambda u: calculate_fresnel_integral(0, u)
    us_to_fresnel = np.vectorize(fresnel_from_zero)

    fresnel_values = us_to_fresnel(us)

    return fresnel_values[0], fresnel_values[1], us


def calculate_fresnel_integral(x0, x1, scaling=1):
    fun_c = lambda x: np.cos(scaling * np.pi * x * x / 2)
    fun_s = lambda x: np.sin(scaling * np.pi * x * x / 2)

    integral_c = lambda x0, x1: integrate.quad(fun_c, x0, x1)
    integral_s = lambda x0, x1: integrate.quad(fun_s, x0, x1)

    integrate_c = np.vectorize(integral_c)
    integrate_s = np.vectorize(integral_s)

    c, c_err = integrate_c(x0, x1)
    s, s_err = integrate_s(x0, x1)

    return c, s


def get_diffraction_pattern(wavelength, slit_width, distance):
    x_lower = -100
    x_higher = 100
    data_points = 1000
    scaling = 2 / (wavelength * distance)

    x = np.linspace(x_lower, x_higher, data_points)

    fresnel_offset_calculation = lambda x: calculate_fresnel_integral(x, x + slit_width, scaling)
    get_wavevector = np.vectorize(fresnel_offset_calculation)

    wavevector_vals = get_wavevector(x)
    return wavevector_vals[0], wavevector_vals[1], x


def get_amplitudes_phases(wavevector_c, wavevector_s, x):


def main():
    # Getting the Cornu spiral figure:
    #    cs,ss,us = get_fresnel_integrals(5,data_points=1000)
    #    plt.plot(cs,ss)
    #    plt.show()

    wv_c, wv_s, x = get_diffraction_pattern(1, 10, 100)

    wv_c_squared = np.power(wv_c, 2)
    wv_s_squared = np.power(wv_s, 2)
    amplitudes = []
    phases = []

    for i in range(0, len(wv_c_squared)):
        amplitudes.append(np.sqrt(wv_c_squared[i] + wv_s_squared[i]))
        phases.append(wv_s[i] / wv_c[i])

    plt.plot(x, amplitudes)
    plt.show()


if __name__ == "__main__":
    main()