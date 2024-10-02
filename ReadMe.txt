# PYTHON ASSIGNMENT 1

This repository contains my python script and the data file.

# Description of the code structure:

Firstly I imported all the usefull modules, in my case matplotlib.pyplot, numpy, and scipy.optimize.

Then I defined the functions.

    The first set of functions is the mathematical functions, that will be used to implement the models. I thus created a second order polynomial function,
    a gaussian function, and a linear combination of the two latter.

    The second set contains six functions. Three of them are diplay functions. The first one returns the parameter best fit values and uncertainties, with 
    the arrays returned by the curve.fit function. The difference with the second one is that it creates an array of values and errors with the optimal values
    and the covariance results of thecurve.fit function, whereas the other one is directly given the list of values and errors. The third one is used to display
    the polynomial function. Then the three remaining functions are defined for three different objectives: first one to find the index in a list of the 
    interval bounds values. For instance in the list [1, 5, 5.4, 8, 99], if I want the index for the interval which bounds are 5 and 8, it will return [1, 3], 
    which are the index of the 5 and 8 numbers. Second one is used to return a given list with a certain part of it masked. That is to say that a certain interval
    of values in the list is going to be replaced by None. Finally the last function is the one that gives initial parameters. This one is usefull for the gaussian
    function, which needs initial parameters values for the curve.fit function to find best estimates.

Then we get to the main function definition.

    First in this function I define a function that is the linear combination function applied with certain parameters fixed, namely the parameters estimations
    for the gaussian and polynomial functions.
    Then I extract the data from the text file and put it in a dictionnary. I create a set of arrays from this dictionnary, that are usefull for certain parts
    of the code. For instance I create arrays which contain data for the peak only, or for the background only. 
    Then I fit the functions and display the results, which consist on three different figures:
        - the full spectrum
        - the full spectrum with the continuum over-plotted
        - the full spectrum with the continuum and the gaussian model over-plotted
    Eventually, I display the results of the parameters obtained with the curve.fit function. 


