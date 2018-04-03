import numpy as np

from astropy.modeling import models, fitting
from matplotlib import pyplot as plt
import matplotlib as mpl

def continuum_normalization(spectrum, absorption=True):
    #-----------------------
    #Plot the spectrum
    #-----------------------
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(spectrum.wave, spectrum.flux)
    ax.set_xlabel('Wavelength $\AA$')
    ax.set_ylabel('Flux')
    
    #-----------------------
    #Normalize the continuum
    #-----------------------
    #Select the edges of the feature and fit a line
    input('Zoom in on line, press enter to continue')
    print('Click on left and right continuum points')
    (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
    x_cont = np.array([x1, x2])
    y_cont = np.array([y1, y2])
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    continuum_fit = np.polyfit(x_cont, y_cont, 1)
    continuum_val = np.polyval(continuum_fit, spectrum.wave)
    #replot the normalized spectrum
    ax.cla()
    flux_norm = spectrum.flux/continuum_val
    ax.plot(spectrum.wave, flux_norm)
    ylim_norm = np.polyval(continuum_fit, xlim)
    ax.set_xlim(x_cont.min()-5, x_cont.max()+5)
    ax.set_ylim(ylim/ylim_norm)
    ax.set_ylim(ymax=1.01)
    #Plot the "continuum" region
    y_cont_norm = y_cont/np.polyval(continuum_fit, x_cont)
    ax.plot(x_cont, y_cont_norm, marker='o', ls='-')
    plt.draw()
    
    #-----------------------
    #Fit the feature
    #-----------------------
    #Select the region to fit
    input('Zoom in on line, press enter to continue') #Add in option to redefine continuum
    print('select left and right edges of fit')
    (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
    #Select the normalized fit spectrum
    line_wave = spectrum.wave[(spectrum.wave<max(x1, x2)) & (spectrum.wave > min(x1, x2))]
    line_flux = flux_norm[(spectrum.wave<max(x1, x2)) & (spectrum.wave > min(x1, x2))]
    #Select the line centers
    print('Select the center(s) of the line(s), right click to continue')
    (x_center, y_center) = plt.ginput(1, timeout=0, show_clicks=True)[0]
    
    #Select the fitting function
    fit_type = input('What shape line would you like to fit? (g=gaussian (default), l=lorentz, m=moffat) ')
    if fit_type == 'l' and absorption is True:
        model = models.Const1D(1.) - models.Lorentz1D(x_0 = x_center)
    elif fit_type == 'm' and absorption is True:
        model = models.Const1D(1.) - models.Moffat1D(x_0 = x_center)
    elif absorption is True:
        model = models.Const1D(1.) - models.Gaussian1D(mean = x_center)
    model.amplitude_0.fixed=True
    
    fitter = fitting.LevMarLSQFitter()
    fit = fitter(model, line_wave, line_flux)
    
    plt.plot(line_wave, line_flux)
    plt.plot(line_wave, fit(line_wave))   
    