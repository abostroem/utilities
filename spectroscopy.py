from astropy import convolution
import numpy as np
import math

def degrade_resolution(flux, disp_in, disp_out):
    '''
    Convolve spectrum with a 1D gaussian kernel to degrade dispersion from disp_in to disp_out
    '''
    fwhm = math.sqrt(disp_out**2 - disp_in**2)
    stddev = fwhm/(2.*math.sqrt(2.*math.log(2.)))
    gauss_kernel=convolution.Gaussian1DKernel(stddev=stddev)
    degraded_flux = convolution.convolve(flux, gauss_kernel)
    return degraded_flux
    
#def test_degrade_resolution