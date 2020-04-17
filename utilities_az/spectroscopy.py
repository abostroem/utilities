from collections import namedtuple

from astropy import convolution
import numpy as np
import math

import extinction

from .define_filters import define_filters, get_cenwave


#Like a light weight class with attributes wavelength and flux
#spectrum1d = namedtuple("spectrum1d", ['wave', 'flux'])
class spectrum1d(object):
    def __init__(self, wave, flux, error=None):
        #Could be np.asarray, but change fit_spectral_lines.py
        self.wave = np.float_(wave)
        self.flux = np.float_(flux)
        self.error = np.float_(error)
        
def degrade_resolution(flux, disp_in, disp_out):
    '''
    Convolve spectrum with a 1D gaussian kernel to degrade dispersion from disp_in to disp_out
    '''
    fwhm = math.sqrt(disp_out**2 - disp_in**2)
    stddev = fwhm/(2.*math.sqrt(2.*math.log(2.)))
    gauss_kernel=convolution.Gaussian1DKernel(stddev=stddev)
    degraded_flux = convolution.convolve(flux, gauss_kernel)
    return degraded_flux
    
def bin_data(data_arr, bin_size=None, bins=[], sum=False, median=False):
    '''
    Bin 1D array of data
    bin_size - this keyword can be set to bin data to a certain number of points. The last point
            may have fewer points
    bins - this keyword can be set to an array of bin edges. Values will always be treated as the left bin edge
    sum - if True, the returned result is the sum of all points in each bin
    median - if True, the returned result is the median of all points in each bin
    
    Note: one and only one of sum and median must be set to True
          one and only one of bin_size and bins must be set
    '''
    binned_array = []
    if bin_size:
        assert len(bins) == 0, 'cannot set both bin_size and bins'
        assert isinstance(bin_size, int), 'bin_size must be an integer'
        arr_len = len(data_arr)
        nbins = int(np.ceil(arr_len/bin_size))
        
        for bin in range(nbins):
            if sum:
                assert median is False, 'cannot set both median and sum to be True'
                binned_array.append(np.sum(data_arr[bin*bin_size:(bin+1)*bin_size]))
            if median:
                assert sum is False, 'cannot set both median and sum to be True'
                binned_array.append(np.median(data_arr[bin*bin_size:(bin+1)*bin_size]))
        return np.array(binned_array)
    if len(bins) != 0:
        assert bin_size is None, 'cannot set both bin_size and bins'
        for left_bin_edge, right_bin_edge in zip(bins[:-1], bins[1:]):
            assert isinstance(left_bin_edge, int) or isinstance(left_bin_edge, np.int64), 'bin edges must be integers'
            if sum:
                assert median is False, 'cannot set both median and sum to be True'
                binned_array.append(np.sum(data_arr[left_bin_edge:right_bin_edge]))
            if median:
                assert sum is False, 'cannot set both median and sum to be True'
                binned_array.append(np.median(data_arr[left_bin_edge:right_bin_edge]))
                
        return np.array(binned_array)
        
def calc_wavelength(header, pixels):
    if ('CTYPE1' in header.keys()) and (len(header['CTYPE1'])>1):
        if header['CTYPE1'] != 'LINEAR': 
            keep_going = input('Attempting to apply a linear solution to non-linear parameters. Continue anyways? (y), n ')
            if keep_going == 'n':
                sys.exit()
    else:
        print('WARNING: No solution type specified, assuming linear')
    CRVAL1 = header['CRVAL1']
    if 'CRPIX1' in header.keys():
        CRPIX1 = header['CRPIX1']
    else:
        CRPIX1 = 0
    if 'CD1_1' in header.keys():
        CD1_1 = header['CD1_1']
    elif 'CDELT1' in header.keys():
        CD1_1 = header['CDELT1']
    elif 'PC1_1' in header.keys():
        CD1_1 = header['PC1_1']
    else:
        sys.exit('Could not identify linear wavelength term (tried CD1_1, CDELT1, PC1_1)')
    
    wavelength = CRVAL1 + CD1_1*(pixels - CRPIX1)
    return wavelength

def apply_redshift(wl, redshift, redden=False):
    '''
    Redden or deredden a wavelength (array) with a given redshift
    
    Inputs:
        wl: wavelength (can be single value or array)
        redshift: redshift value
        redden: If False, code will calculate rest wavelength assuming input is observed
                If True, code will calculate observed wavelength assuming input is rest
    '''
    if isinstance(wl, list):
        wl = np.array(wl)
    if redden is False:
        new_wl = wl/(1+redshift)
    elif redden is True:
        new_wl = wl * (1+redshift)
    return new_wl
    
def scale_spectra(spec1, spec2, wlmin = None, wlmax=None, scale_factor=False):
    '''
    Scale spec1 to spec2 by integrating from wlmin to wlmax and using
    the area as a scaling factor
    Inputs:
        spec1, spec2: two spectrum1d objects with attributes wave and flux
        wlmin: shortest wavelength to be used in scaling
        wlmax: longest wavelength to be used in scaling
    Outputs:
        spec1: spec1d object with scaled flux
        scale_factor (optional): if scale_factor = True, then the multiplicative factor is returned
    '''    
    if wlmin is None:
        wlmin = max(spec1.wave.min(), spec2.wave.min())
    if wlmax is None:
        wlmax = min(spec1.wave.max(), spec2.wave.max())
    windx_1 = (spec1.wave >= wlmin) &  (spec1.wave <= wlmax)
    windx_2 = (spec2.wave >= wlmin) & (spec2.wave <= wlmax)
    delta_wl1 = np.mean(spec1.wave[1:] - spec1.wave[:-1])
    delta_wl2 = np.mean(spec2.wave[1:] - spec2.wave[:-1])
    area1 = np.sum(spec1.flux[windx_1])*delta_wl1
    area2 = np.sum(spec2.flux[windx_2])*delta_wl2
    spec1 = spectrum1d(spec1.wave, spec1.flux/area1*area2)
    if scale_factor is True:
        return spec1, area2/area1
    else:
        return spec1
    
def correct_for_galactic_extinction(spec, E_BV, R_V=3.1):
    '''
    Correct flux for galactic (Milky Way) extinction using a Cardelli law
    Inputs:
        spec: spectrum1d object to be corrected (has wave and flux attributes)
        E_BV: E(B-V) for correction
        R_V: default to 3.1 values
    Output:
        spectrum1d object with the dust corrected flux
    '''
    A_V = R_V*E_BV
    new_flux = extinction.apply(-extinction.ccm89(spec.wave, A_V, R_V), spec.flux)    
    return spectrum1d(spec.wave, new_flux)   
    
def calc_extinction(E_BV, band, Rv=3.1, E_BV_err=0):
    '''
    Calculate extinction using the Cardelli law
    E_BV: float
        E(B-V) value
    band: str
        filter corresponding to filter in define_filters dictionary
    Rv: float
        R(V) set to 3.1 by default
    E_BV_err: float
        1 sigma error in E(B-V), set to 0 by default
    '''
    Av = Rv*E_BV
    Av_err = Rv**2 * E_BV_err**2
    cenwave = np.array([float(get_cenwave(band))])
    A_band = extinction.ccm89(cenwave, Av, Rv)
    if Av == 0:
        A_band_err = 0.0
    else:
        A_band_err = np.sqrt(Av_err * A_band/Av)
    return A_band, A_band_err

############################
#Tests
############################
def test_full_bins():                
    test_arr = np.arange(100)
    binned_array1 = bin_data(test_arr, bin_size = 10, sum=True)
    assert (binned_array1 == np.arange(45, 1000, 100)).all(), 'failed test on full bins'

def test_incomplete_last_bin():
    test_arr = np.arange(100)
    binned_array2 = bin_data(test_arr, bin_size = 15, sum=True)
    assert (binned_array2 == np.array([105, 330, 555, 780, 1005, 1230, 945])).all(), 'failed test on last bin not being full'

def test_premade_bins():
    test_arr = np.arange(100)
    bins = np.arange(0, 101, 10)
    binned_array1 = bin_data(test_arr, bin_size = 10)
    binned_array3 = bin_data(test_arr, bins = bins)
    assert (binned_array3 == binned_array1).all(), 'bins not functioning as intended'
    
def test_scale_spectra():
    spec1 = spec.spectrum1d(np.arange(6000, 8000, 1), np.ones(2000))
    spec2 = spec.spectrum1d(np.arange(6000, 8000, 1), np.ones(2000)*4.2)
    new_spec, scale_factor = spec.scale_spectra(spec1, spec2, scale_factor=True)
    plt.figure()
    plt.plot(spec1.wave, spec1.flux, ':', label='original flux')
    plt.plot(spec1.wave, spec1.flux*scale_factor, label='scaled flux')
    plt.plot(spec2.wave, spec2.flux, '--', label='flux to scale to')
    plt.legend(loc='best')
    assert np.isclose(scale_factor,4.2), 'scale_factor={}'.format(scale_factor)
    
def test_scale_spectra2():
    flux1 = np.ones(2000)
    flux1[1000:1500]+=1
    flux2 = np.ones(2000)+4.2
    flux2[1000:1500]+=3
    spec1 = spec.spectrum1d(np.arange(6000, 8000, 1), flux1)
    spec2 = spec.spectrum1d(np.arange(6000, 8000, 1), flux2)
    new_spec, scale_factor = spec.scale_spectra(spec1, spec2, scale_factor=True, wlmin=7000, wlmax=7500)
    plt.figure()
    plt.plot(spec1.wave, spec1.flux, ':', label='original flux')
    plt.plot(spec1.wave, spec1.flux*scale_factor, label='scaled flux')
    plt.plot(spec2.wave, spec2.flux, '--', label='flux to scale to')
    plt.legend(loc='best')
    assert np.isclose(scale_factor,4.101098901098902), 'scale_factor={}'.format(scale_factor)


def test_scale_spectra3():
    flux1 = np.ones(2000)
    flux1[1000:1500]+=1
    flux2 = np.ones(2000)+4.2
    flux2[1000:1500]+=3
    spec1 = spec.spectrum1d(np.arange(6000, 8000, 1), flux1)
    spec2 = spec.spectrum1d(np.arange(6000, 8000, 1), flux2)
    new_spec, scale_factor = spec.scale_spectra(spec1, spec2, scale_factor=True, wlmin=7200, wlmax=7250)
    plt.figure()
    plt.plot(spec1.wave, spec1.flux, ':', label='original flux')
    plt.plot(spec1.wave, spec1.flux*scale_factor, label='scaled flux')
    plt.plot(spec2.wave, spec2.flux, '--', label='flux to scale to')
    plt.legend(loc='best')
    assert np.isclose(scale_factor,4.101098901098902, rtol=0.1), 'scale_factor={}'.format(scale_factor)
    


