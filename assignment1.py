import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# ----------------------------- MODELS FUNCTIONS DEFINITION ------------------------

def gaussian_model(x, mu, sig, A, c):
    """
    returns a gaussian function (amplitude non-normalized)
    """
    return A * np.exp( - (x - mu) ** 2 / (2 * sig ** 2) ) + c

def polynomial_model(x, a, b, c):
    """
    returns a polynomial function
    """
    return a * x ** 2 + b * x + c

def combination_model(x, mu, sig, A, d, a, b, c, alpha, beta):
    """ 
    returns a linear combination of the gaussian and polynomial functions
    """
    return alpha * gaussian_model(x, mu, sig, A, d) + beta * polynomial_model(x, a, b, c)


# ---------------------------------- OTHER FUNCTIONS ------------------------

def format_result(params, popt, pcov):
    """
    Displays parameter best estimates and uncertainties
    """
    perr = np.sqrt(np.diag(pcov))
     
    _lines = (f"{p} = {round(o, 4)} ± {round(e, 4)}" for p, o, e in zip(params, popt, perr))
    return "\n".join(_lines)

def format_simple(params, value, error):    #specific for the compound function, when I have to give values and errors directly
    """                                     
    Display parameter best estimates and uncertainties when values are given
    """
    _lines = (f"{p} = {round(o, 4)} ± {round(e, 4)}" for p, o, e in zip(params, value, error))
    return "\n".join(_lines)

def polynomial_function_display(popt):
    """ 
    Displays parameter best estimates in the polynomial equation format
    """
    return f"{round(popt[0], 4)}x² + {round(popt[1], 2)}x + {round(popt[2], 1)}"


def find_wavelength_peak(list): #knowing the values of the bounds of the interval, gives the interval in terms of index
    ''' 
    gives the index range (in the list) of the peak
    '''
    i = 0
    while list[i] < 6680:       # 6680 and 6695 values determinated with a visual inspection of the spectra
        i += 1
    min = i
    while list[i] < 6695:
        i += 1
    maxx = i
    return [min, maxx]

def only_continuum(list, min_ ,max_):   #used to have only the background part of the raw data, the rest is masked
    """ 
    returns the same list but with a masked interval
    """
    new_list = np.asarray([i for i in list])
    for i in range(min_, max_):
        new_list[i] = None
    return new_list

def initial_parameters(list):   #here we set the amplitude A as the maximum value of the peak, the mu value corresponding to the point on which
                                #the peak is centered is set as the interval center, and sigma is the FWHM of the peak (up to a factor) set to 10
    """
    returns guessed initial parameters for the gaussian function
    """
    maxx = 0
    for i in list:
        if i > maxx:
            maxx = i
    A = maxx
    mu = (6680 + 6695)/2
    sig = 10 / 2 * np.sqrt(2* np.log(2))
    return [mu, sig, A, 0]

# ---------------------------------- MAIN FUNCTION ----------------------------

def main():
    """ 
    returns 3 figures: 
    - the spectrum of the raw data
    - the spectrum with the background fit
    - the spectrum with the compound model (background and peak)
    and the values of the parameters of the fitting function with uncertainties
    """

    def combination_model_applied(x, alpha, beta):  #for the fitting of the linear combination, only changing parameters are alpha and beta
        return combination_model(x, *popt, *popt_poly, alpha, beta)

    filename = 'spectrum.txt'

    data = {'Wavelength' : [], 'Flux' : []}

    with open(filename, 'r') as file:
        lines = file.readlines()[27:]          
        for line in lines:
            wavelength, flux = line.split(',')  
            data['Wavelength'].append(float(wavelength))    
            data['Flux'].append(float(flux))
    
    # data arrays:
    wavelength = np.asarray(data['Wavelength'])
    flux = np.asarray(data['Flux'])
    peak = find_wavelength_peak(wavelength)         # index interval of the peak
    flux_continuum = only_continuum(flux, *peak)    # flux data without the peak part
    wavelength_peak = wavelength[peak[0] : peak[1]] # wavelenghts of the peak only
    flux_peak = flux[peak[0] : peak[1]]             # fluxes of the peak only

    initial_param = initial_parameters(flux_peak)

    #fitting the functions (first the background with the polynomial one, then the peak with the gaussian one, and eventually the linear combination )
    popt_poly, pcov_poly = curve_fit(polynomial_model, wavelength, flux_continuum, nan_policy='omit')   # "nan_policy='omit'" means that the calculation 
    np.sqrt(np.diag(pcov_poly))                                                                         # are performed ignoring the NaN values
    popt, pcov = curve_fit(gaussian_model, wavelength_peak, flux_peak, p0=initial_param)
    np.sqrt(np.diag(pcov))
    popt_combination, pcov_combination = curve_fit(combination_model_applied, wavelength, flux)
    np.sqrt(np.diag(pcov_combination))

    #plot the results:
    # -- the full spectrum
    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(wavelength, flux, color = 'darkslateblue', label = 'raw data')

    plt.title('Spectrum', size=17)
    plt.xlabel('wavelength (Å)')
    plt.ylabel('flux (ADU)')
    plt.legend()

    # -- the full spectrum with the continuum over-plotted
    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(wavelength, flux, color = 'darkolivegreen', label = 'raw data')
    poly_func = polynomial_function_display(popt_poly)
    ax.plot(wavelength, polynomial_model(wavelength, *popt_poly), c = 'tomato', label = f'polynomial model:\n f(x) = {poly_func}')

    plt.title('Continuum of the spectrum', size=17)
    plt.xlabel('wavelength (Å)')
    plt.ylabel('flux (ADU)')
    plt.legend()

    # -- the full spectrum with the compound model over-plotted
    fig, ax = plt.subplots(figsize = (10, 5))
    ax.plot(wavelength, flux, color = 'cadetblue', label = 'raw data')
    ax.plot(wavelength, combination_model_applied(wavelength, *popt_combination), c = 'firebrick', label = 'compound model')

    # in order to display the values of the final parameters, some calculations are needed (for instance new_amplitude = old_one * alpha):
    # we also need to multiply some of the uncertainties to have the final one for some parameters such as the amplitude
    values = (popt[0], popt[1], popt[2]*popt_combination[0], popt_combination[1]*popt_poly[0], popt_combination[1]*popt_poly[1],\
            popt_combination[1]*popt_poly[2]+popt[3]*popt_combination[0])
    
    perr_gaussian = np.sqrt(np.diag(pcov))    
    perr_poly = np.sqrt(np.diag(pcov_poly))    
    perr_compound = np.sqrt(np.diag(pcov_combination))   

    errors = np.array([perr_gaussian[0], perr_gaussian[1], perr_gaussian[2]*perr_compound[0], perr_poly[0]*perr_compound[1], \
            perr_poly[1]*perr_compound[1], perr_poly[2]*perr_compound[1]])
    
    # expression of the final function
    model = f'{round(values[0], 3)} * exp(-(x-{round(values[1], 3)})²/2*{round(values[2], 3)}) + {round(values[3], 3)}x² +{round(values[4], 3)}x + {round(values[5], 3)}'

    plt.title('Spectrum compound fitting', size=17)
    plt.xlabel('wavelength (Å)')
    plt.ylabel('flux (ADU)')
    plt.text(6640, 7, f"Compound function : f(x) = {model}", fontsize=8, color='firebrick')  #display the equation of the model on the figure

    plt.legend()
    plt.show()

    # display all the parameters values. First the values for the background fitting only, then for the gaussian fitting of the peak only, 
    # and eventually the linear combination fitting and the whole model
    params_poly = ('a', 'b', 'c')
    params_gaussian = ('µ', 'σ', 'A', 'c')
    params_compound = ('α', 'β')
    params = ('µ', 'σ', 'A', 'a', 'b', 'c')

    print(' parameters values, for the fitting function f(x) = ax² + bx + c: \n')
    print(format_result(params_poly, popt_poly, pcov_poly))
    print('\n parameters values, for the fitting function f(x) = A * exp(-(x-µ)²/2σ²) + d \n')
    print(format_result(params_gaussian, popt, pcov))
    print('\n parameters values, for the linear combination f(x) = α * gaussian_model + β * polynomial_model \n')
    print(format_result(params_compound, popt_combination, pcov_combination))
    print('\n parameter values for the fitting function f(x) = A * exp(-(x-µ)²/2σ²) +  ax² + bx + c : \n')
    print(format_simple(params, values, errors))

if __name__ == "__main__":

    text = (
    
        """\n * * * Method: * * *\n The first graph represents the data from the text document directly. We can clearly see that there is two part in this spectrum:"""
        """There is the background and the peak, in the [6680 : 6695] interval. The objective is to create a model of the spectrum, and for that, we """
        """are using the fact that this spectrum can be broken down in two parts. First we fit the background, then the peak, and eventually a linear """
        """combination of the two.\n\n - to fit the background, I used a second order polynomial function. The results shows a very flat parabola, so a """
        """linear function would have certainly worked too. \n - to fit the peak I used a gaussian function, initializing the parameters with the amplitude""" 
        """A being the maximum of the peak, the centre of the gaussian function µ the centre of the interval I defined as the peak interval, and σ which """
        """is up to a factor (2 sqrt(2 ln(2))) the FWHM of the peak at 10 (*2 sqrt(2 ln(2))).\n - for the combined model, I fitted a linear combination of"""
        """the two previous fitting functions (with the parameters set at the values of the best fit)\n * * * \n\n"""
    )
    print(text)
    main()