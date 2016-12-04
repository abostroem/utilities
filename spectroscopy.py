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
    assert header['CTYPE1'] == 'LINEAR'
    CRVAL1 = header['CRVAL1']
    CRPIX1 = header['CRPIX1']
    CD1_1 = header['CD1_1']
    
    wavelength = CRVAL1 + CD1_1*(pixels - CRPIX1)
    return wavelength


############################
#Tests
############################
def test_full_bins():                
    test_arr = np.arange(100)
    binned_array1 = bin_data(test_arr, bin_size = 10)
    assert (binned_array1 == np.arange(45, 1000, 100)).all(), 'failed test on full bins'

def test_incomplete_last_bin():
    binned_array1 = bin_data(test_arr, bin_size = 10)
    binned_array2 = bin_data(test_arr, bin_size = 15)
    assert (binned_array2 == np.array([105, 330, 555, 780, 1005, 1230, 945])).all(), 'failed test on last bin not being full'
def test_premade_bins():
    bins = np.arange(0, 101, 10)
    binned_array1 = bin_data(test_arr, bin_size = 10)
    binned_array3 = bin_data(test_arr, bins = bins)
    assert (binned_array3 == binned_array1).all(), 'bins not functioning as intended'