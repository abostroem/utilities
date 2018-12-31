import os
import shutil
import sys

#Anaconda
import numpy as np
from astropy.time import Time, TimeDelta
from astropy import units as U
from astropy.io import ascii as asc
from astropy.io import fits

#Pip
import extinction
import pymysql
#Local
from .define_filters import define_filters
from . import connect_to_sndavis


class LightCurve(object):
    def __init__(self, name, ra, dec, type=None, dist_mod=None, dist_mod_err=None, 
                       ebv_mw=None, ebv_mw_err=None, ebv_host=None, ebv_host_err=None, 
                       jdexpl=None, jdexplerr=None):
        self.name=name
        self.ra = result['ra0']
        self.dec = result['dec0']
        self.type = result['sntype'] #TODO - update this to look at the type table and make a human readable type
        self.dist_mod = result['mu']
        self.dist_mod_err = result['muerr']
        self.ebv_mw = result['ebvg']
        self.ebv_mw_err = result['ebvgerr']
        self.ebv_host = result['ebvi'] 
        self.ebv_host_err = result['ebvierr']
        self.jdexpl = result['jdexpl']
        self.jdexpl_err = result['jdexplerr']
        self.phase = {}
        self.jd = {}
        self.apparent_mag = {}
        self.apparent_mag_err = {}
        self.abs_mag = {}
        self.abs_mag_err = {}
    def add_photometry(filter, jd, phot, phot_err):
        self.jd[filter] = jd
        self.apparent_mag[filter] = phot
        self.apparent_mag_err[filter] = phot_err
        if self.jdexpl:
            self.phase[filter] = self.jd[filter] -self.jdexpl
        
        
        
def get_cenwave(band):
    band_dict = define_filters()
    if band in band_dict.keys():
        band_info = band_dict[band]
        cenwave = band_info[2]
        return cenwave
    else:
        print('band {} not in band dict keys: {}'.format(band, band_dict.keys()))
        raise 'KeyError' #TODO: figure out why this doesn't work

def calc_abs_mag(app_mag, dist_mod, mw_ebv, band, host_ebv=0, dist_mod_err=0, app_mag_err=0):
    R_v = 3.1
    A_v_mw = R_v*mw_ebv
    A_v_host = R_v* host_ebv
    cenwave = np.array([float(get_cenwave(band))])
    A_band_mw = extinction.ccm89(cenwave, A_v_mw, R_v)
    A_band_host = extinction.ccm89(cenwave, A_v_host, R_v)
    abs_mag = app_mag - dist_mod - A_band_mw - A_band_host
    abs_mag_err = np.sqrt(app_mag_err**2 + dist_mod_err**2) #check reddening error
    return abs_mag, abs_mag_err

class LightCurve2(object):
    '''
    Retrieve information from sndavis database
    '''
    def __init__(self, name):
        self.name=name
        self.db, self.cursor = connect_to_sndavis.get_cursor()
        self.id = self.get_sn_id()
        if self.id is not None:
            self.cursor.execute('SELECT * FROM idsupernovae WHERE id = {}'.format(self.id))
            result = self.cursor.fetchone()
            self.ra = result['ra0']
            self.dec = result['dec0']
            self.type = result['sntype'] #TODO - update this to look at the type table and make a human readable type
            self.dist_mod = result['mu']
            self.dist_mod_err = result['muerr']
            self.ebv_mw = result['ebvg']
            self.ebv_mw_err = result['ebvgerr']
            self.ebv_host = result['ebvi'] 
            self.ebv_host_err = result['ebvierr']
            self.quality = result['quality']
            self.jdexpl = result['jdexpl']
            self.jdexpl_err = result['jdexplerr']
            self.jd = {}
            self.apparent_mag = {}
            self.apparent_mag_err = {}
            self.phase = {}
            self.abs_mag = {}
            self.abs_mag_err = {}
            self.slopes = {'s1':{}, 's2':{}, 's50':{}, 'tail':{}, 
                            's1_range':{}, 's2_range':{}, 's50_range':{}, 'tail_range':{},
                            's1_err':{}, 's2_err':{}, 's50_err':{}, 'tail_err':{}}
        else:
            print('SN {} not found in database'.format(self.name))
            #TODO add an absolute magnitude error
             
    def get_sn_id(self):
        self.cursor.execute('SELECT targetid FROM supernovanames WHERE name = {}'.format('"{}"'.format(self.name)))
        result = self.cursor.fetchone()
        if result is None:
            return None
        return result['targetid']
 
    def get_photometry(self, band='all'):
        '''
        Get apparent magnitude of a specific band
        '''
        if band == 'all':
            self.cursor.execute("SELECT DISTINCT BINARY(filter) FROM photometry WHERE targetid={}".format(self.id))
            results = self.cursor.fetchall()
            bands = [iband['BINARY(filter)'].decode('utf-8') for iband in results]

        elif len(band)==1: #Single filter
            bands = [band]
        else: #multiple filters
            bands = band
        for ifilter in bands:   
            self.cursor.execute("SELECT jd, mag, magerr FROM photometry WHERE targetid={} AND filter=BINARY('{}')".format(self.id, ifilter))
            results = self.cursor.fetchall()
            jd = []
            mag = []
            mag_err = []
            for irow in results:
                jd.append(float(irow['jd']))
                mag.append(float(irow['mag']))
                mag_err.append(float(irow['magerr']))
            self.jd[ifilter] = np.array(jd)
            self.apparent_mag[ifilter] = np.array(mag)
            self.apparent_mag_err[ifilter] = np.array(mag_err)
            self.phase[ifilter] = self.jd[ifilter] - self.jdexpl
            
    def get_abs_mag(self, band='all'):
        '''
        Calcualte absolute magnitude
        '''
        if band == 'all':
            bands = self.apparent_mag.keys()
        elif len(band)==1:
            bands = [band]
        else:
            bands = band
        if self.dist_mod_err is None:
            self.dist_mod_err = 0
        for iband in bands:
            if self.ebv_host is None:
                self.abs_mag[iband], self.abs_mag_err[iband] = calc_abs_mag(self.apparent_mag[iband], 
                                                   self.dist_mod, 
                                                   self.ebv_mw, iband, 
                                                   app_mag_err=self.apparent_mag_err[iband], 
                                                   dist_mod_err=self.dist_mod_err)
            else:
                self.abs_mag[iband], self.abs_mag_err[iband] = calc_abs_mag(self.apparent_mag[iband], 
                                                                            self.dist_mod, 
                                                                            self.ebv_mw, 
                                                                            iband, 
                                                                            host_ebv=self.ebv_host,
                                                                            app_mag_err=self.apparent_mag_err[iband], 
                                                                            dist_mod_err=self.dist_mod_err)
    
    def get_slope(self,slope_type, band='V'):
        '''
        '''
        self.cursor.execute('SELECT * FROM snslope WHERE  targetid={} AND \
                                                         slopetype="{}" AND \
                                                         filter=BINARY("{}")'.format(self.id, slope_type, band))
        result = self.cursor.fetchone()
        if result is not None:
            self.slopes['{}'.format(slope_type)][band]       = result['slope']
            self.slopes['{}_err'.format(slope_type)][band]   = result['slopeerr']
            self.slopes['{}_range'.format(slope_type)][band] = (result['tstart'], result['tstop'])
        else:
            self.slopes['{}'.format(slope_type)][band]       = None
            self.slopes['{}_err'.format(slope_type)][band]   = None
            self.slopes['{}_range'.format(slope_type)][band] = [None, None]


def get_closest_photometry(date_obs, jd, phot_mag):
    '''
    Find the closest photometric observations to a given date
    
    Inputs:
    ----------
        date_obs: Time object
            astropy Time object of the observervation date that you want to find the closest observations to
        jd: list
            list of photometric observations in JD
        phot_mag: list
            list of magnitudes corresponding to the jd list for a single filter
        
    Returns:
    -----------
        closest_mag: float
            the closest magnitude to the date_obs given
        date_sep: float
            the difference in days between the observation and the 
    '''
    date_indx = np.argmin(np.abs(jd - date_obs.jd))
    closest_mag = phot_mag[date_indx]
    date_sep = (jd - date_obs.jd)[date_indx]
    return closest_mag, date_sep

def convert_mag_to_flux(phot_mag, ifilter):
    '''
    Converts magnitude units to flux units (ergs/cm^2/s/A)
    
    Inputs:
    --------
    phot_mag: list or array
        list or array of magnitudes
    ifilter: str
        filter that phot_mag observations were taken with
        
    Returns:
    --------
    flux : array
        list of flux values
    '''
    bandpar_dict = define_filters()
    flux = 10**(np.array(phot_mag)/-2.5)*bandpar_dict[ifilter][4]
    return flux
    
def get_cenwave(ifilter):
    '''
    return the central wavelength of a given filter
    
    Inputs:
    -------
    ifilter: str
        the name of a filter (must be present in define_filters)
    
    Returns:
    --------
    cenwave: float
        central wavelength of filter
    '''
    bandpar_dict = define_filters()
    cenwave = bandpar_dict[ifilter][2]
    return cenwave

def scale_spectra_quba(snname, filename, filter_dir=None, date_kw='date-obs', max_sep=7, 
                        header_date=False, sndavis=True, lightcurve=None):
    '''
    You must be connected to dark
    filter_dir should be /Users/bostroem/Dropbox/DLT40_all/script/scalespectra/filters/
    if header_date is True, then date_kw contains the keyword to read from the header.
    if header_date is False, then the date_kw should be treated as the date of the observation
    
    Inputs:
    ---------
    snname: str
        name of supernova - will be used to find in SNDAVIS database is sndavis is True
    filename: str
        name of spectrum file. 
    filter_dir: str
        location of information filter throughput
    date_kw: str
        if header_date is True: name of keyword in fits file to use to find date
        if header_date is False: date of observation
    max_set: float
        number of days that can separate photometry point and date of spectrum
    header_date: bool
        controls whether spectrum date is read from header
    sndavis: bool
        if True: use SNDAVIS database to get light curve
        if False: pass a light curve into the light curve keyword
    lightcurve: LightCurve object
        an object that has a jd attribute and an app_mag attribute that are dictionaries with filters as keys
    
    
    
    '''
    import qubascalespectra
    if sndavis:
        lightcurve = LightCurve2(snname)
        lightcurve.get_photometry()
    if header_date is True:
        date_obs = Time(fits.getval(filename, date_kw, 0))
    else:
        date_obs = Time(date_kw)
    spec_phase = (date_obs - Time(lightcurve.jdexpl, format='jd')).value
    band = ''
    mphot = []
    for ifilter in lightcurve.jd.keys():
        if len(ifilter) == 1: #avoid swift filters
            cenwave = get_cenwave(ifilter)
            if (cenwave > 3500):
                mag, date_sep = get_closest_photometry(date_obs, lightcurve.jd[ifilter], lightcurve.apparent_mag[ifilter])
                if (date_sep <= max_sep):
                    band=band+ifilter
                    mphot.append(mag)
    if len(mphot) > 1:
        qubascalespectra.scale_spectrum(filename, band, mphot, filter_dir)
    else:
        print('WARNING: no photometry within {} days of observation date ({})'.format(max_sep, date_obs))
    shutil.move('log.txt', filename.split('.')[0]+'.log')
    with open(filename.split('.')[0]+'.log', 'a') as ofile:
        ofile.write(band+'\n')
        mphot_str = [str(iphot) for iphot in mphot]
        ofile.write(','.join(mphot_str)+'\n')
        
            
    
        
