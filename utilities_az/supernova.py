import os
import shutil
import sys
from collections import defaultdict
import warnings

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
from .define_filters import define_filters, get_cenwave
from . import connect_to_sndavis, connect_to_supernova
from . import spectroscopy as spec


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
        self.phase = defaultdict(list)
        self.jd = defaultdict(list)
        self.apparent_mag = defaultdict(list)
        self.apparent_mag_err = defaultdict(list)
        self.abs_mag = defaultdict(list)
        self.abs_mag_err = defaultdict(list)
        self.A_host = defaultdict(list)
        self.A_err_host = defaultdict(list)
        self.A_mw = defaultdict(list)
        self.A_err_mw = defaultdict(list)
        
    def add_photometry(filter, jd, phot, phot_err):
        self.jd[filter] = jd
        self.apparent_mag[filter] = phot
        self.apparent_mag_err[filter] = phot_err
        if self.jdexpl:
            self.phase[filter] = self.jd[filter] -self.jdexpl

def calc_abs_mag(app_mag, dist_mod, A_mw, A_host=0, dist_mod_err=0, app_mag_err=0, A_err_mw=0, A_err_host=0):
    '''
    Calculate the absolute magnitude in a given filter given the apparent magnitude, distance modulus, and extinction
    app_mag: arr or float
        array of apparent magnitudes
    dist_mod: float
        distance modulus
    A_mw: float
        Magnitudes of extinction in filter due to the Milky Way
    A_host: float
        Magnitudes of extinction in filter due to host; default=0
    dist_mod_err: float
        error in the distance modulus; default=0
    app_mag_err: arr or float
        array of error values for apparent magnitude (should be same size as app_mag); default=0
    A_err_mw: float
        error on A_mw (in magnitudes); default=0
    A_err_host: float
        err on A_host (in magnitudes); default=0
    '''
    abs_mag = app_mag - dist_mod - A_mw - A_host
    abs_mag_err = np.sqrt(app_mag_err**2 + dist_mod_err**2 + A_err_mw**2 + A_err_host**2)
    return abs_mag, abs_mag_err

class LightCurve2(object):
    '''
    Retrieve information from sndavis database
    
    Absolute magnitudes are corrected for extinction, apparent magnitude is not
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
            if result['ebvi'] is None:
                self.ebv_host = 0
                self.ebv_host_err = 0
                print('WARNING: host extinciton is None, setting to 0')
            else:
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
            self.A_host = {}
            self.A_err_host = {}
            self.A_mw = {}
            self.A_err_mw = {}
            self.jd_limit = {}
            self.phase_limit = {}
            self.apparent_mag_limit = {}
            self.apparent_mag_err_limit = {}
            self.abs_mag_limit = {}
            self.abs_mag_err_limit = {}
            self.telescope = {}
            self.telescope_limit = {}
            self.project = {}
            self.project_limit = {}
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
        filter_dict =  define_filters()
        if band == 'all':
            self.cursor.execute("SELECT DISTINCT BINARY(filter) FROM photometry WHERE targetid={}".format(self.id))
            results = self.cursor.fetchall()
            bands = [iband['BINARY(filter)'].decode('utf-8') for iband in results]

        elif len(band)==1 or isinstance(band, str): #Single filter
            bands = [band]
        else: #multiple filters
            bands = band
        for ifilter in bands:   
            if ifilter in filter_dict.keys():
                A_host_band, A_err_host_band = spec.calc_extinction(self.ebv_host, ifilter)
                A_mw_band, A_err_mw_band = spec.calc_extinction(self.ebv_mw, ifilter)
            else:
                print('Could not calculate extinction for {} because not in filter dictionary'.format(ifilter))
                A_host_band = 0
                A_err_host_band = 0
                A_mw_band = 0
                A_err_mw_band = 0
            #self.cursor.execute("SELECT jd, mag, magerr, datatype, configuration  FROM photometry JOIN photometrysource on photometry.source=photometrysource.source WHERE targetid={} AND filter=BINARY('{}')".format(self.id, ifilter))
            self.cursor.execute("SELECT jd, mag, magerr, datatype, source FROM photometry WHERE targetid={} AND filter=BINARY('{}')".format(self.id, ifilter))
            results = self.cursor.fetchall()
            jd = []
            mag = []
            mag_err = []
            telescope = []
            project=[]
            jd_limit = []
            mag_limit = []
            mag_err_limit = []
            telescope_limit = []
            project_limit = []
            for irow in results:
                if irow['magerr'] is None:
                    irow['magerr'] = 0
                if irow['datatype'] < 0: #photometry represents a limit
                    jd_limit.append(float(irow['jd']))
                    mag_limit.append(float(irow['mag']))
                    mag_err_limit.append(float(irow['magerr']))
                    #telescope_limit.append(irow['configuration'])
                    #telescope_limit.append(irow['source'])
                    #if irow['configuration'].startswith('1m'):
                    #    project_limit.append('LCO')
                    #elif irow['configuration'] == 'Meckering' or irow['configuration'] == 'Prompt5':
                    #    project_limit.append('DLT40')
                    #else:
                    #    project_limit.append(None)
                    telescope_limit.append(irow['source'])
                    if irow['source'] in [5209, 5109, 5009, 4809, 4909, 4709, 4609, 4509, 140, 141, 142, 165, 166, 167]:
                        project_limit.append('LCO')
                    elif irow['source'] in [5300, 5400]:
                        project_limit.append('DLT40')
                    else:
                        project_limit.append(None)
                else:
                    jd.append(float(irow['jd']))
                    mag.append(float(irow['mag']))
                    mag_err.append(float(irow['magerr']))
                    # telescope.append(irow['configuration'])
                    # if irow['configuration'].startswith('1m'):
                    #     project.append('LCO')
                    # elif irow['configuration'] == 'Meckering' or irow['configuration'] == 'Prompt5':
                    #     project.append('DLT40')
                    # else:
                    #     project.append(None)
                    telescope.append(irow['source'])
                    if irow['source'] in [5209, 5109, 5009, 4809, 4909, 4709, 4609, 4509, 140, 141, 142, 165, 166, 167]:
                        project.append('LCO')
                    elif irow['source'] in [5300, 5400]:
                        project.append('DLT40')
                    else:
                        project.append(None)
            self.jd[ifilter] = np.array(jd)
            self.apparent_mag[ifilter] = np.array(mag)
            self.apparent_mag_err[ifilter] = np.array(mag_err)
            if self.jdexpl is not None:
                self.phase[ifilter] = self.jd[ifilter] - self.jdexpl
            else:
                print('WARNING: No explosion epoch in database, phase is not calculated')
            self.A_host[ifilter] = A_host_band
            self.A_err_host[ifilter] = A_err_host_band
            self.A_mw[ifilter] = A_mw_band
            self.A_err_mw[ifilter] = A_err_mw_band
            self.jd_limit[ifilter] = np.array(jd_limit)
            self.apparent_mag_limit[ifilter] = np.array(mag_limit)
            self.apparent_mag_err_limit[ifilter] = np.array(mag_err_limit)
            if self.jdexpl is not None:
                self.phase_limit[ifilter] = self.jd_limit[ifilter] - self.jdexpl
            else:
                print('WARNING: No explosion epoch in database, phase is not calculated')
            self.telescope[ifilter] = np.array(telescope)
            self.telescope_limit[ifilter] = np.array(telescope_limit)
            self.project[ifilter] = np.array(project)
            self.project_limit[ifilter] = np.array(project_limit)
            

            
    def get_abs_mag(self, band='all'):
        '''
        Calcualte absolute magnitude
        '''
        print('Calculating Absolute Magntidue with Extinction')
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
                                                   self.A_mw[iband], 
                                                   app_mag_err=self.apparent_mag_err[iband], 
                                                   dist_mod_err=self.dist_mod_err)
                self.abs_mag_limit[iband], self.abs_mag_err_limit[iband] = calc_abs_mag(self.apparent_mag_limit[iband], 
                                                   self.dist_mod, 
                                                   self.A_mw[iband], 
                                                   app_mag_err=self.apparent_mag_err_limit[iband], 
                                                   dist_mod_err=self.dist_mod_err)
            else:
                self.abs_mag[iband], self.abs_mag_err[iband] = calc_abs_mag(self.apparent_mag[iband], 
                                                                            self.dist_mod, 
                                                                            self.A_mw[iband], 
                                                                            A_host=self.A_host[iband],
                                                                            app_mag_err=self.apparent_mag_err[iband], 
                                                                            dist_mod_err=self.dist_mod_err,
                                                                            A_err_mw=self.A_err_mw[iband],
                                                                            A_err_host=self.A_err_host[iband])
                self.abs_mag_limit[iband], self.abs_mag_err_limit[iband] = calc_abs_mag(self.apparent_mag_limit[iband], 
                                                                            self.dist_mod, 
                                                                            self.A_mw[iband], 
                                                                            A_host=self.A_host[iband],
                                                                            app_mag_err=self.apparent_mag_err_limit[iband], 
                                                                            dist_mod_err=self.dist_mod_err,
                                                                            A_err_mw=self.A_err_mw[iband],
                                                                            A_err_host=self.A_err_host[iband])
    
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

########################################

class LightCurveLCO(object):
    '''
    Retrieve information from LCO supernova database
    '''
    def __init__(self, name, **kwargs):
        if 'config_dir' not in kwargs.keys():
            kwargs['config_dir'] = os.environ['LCOSNDIR']
        if 'config_filename' not in kwargs.keys():
            kwargs['config_filename'] = 'configure'
        self.name=name
        self.db, self.cursor = connect_to_supernova.get_cursor(**kwargs)
        self.id = self.get_sn_id()
        if self.id is not None:
            self.cursor.execute('SELECT * FROM targets WHERE id = {}'.format(self.id))
            result = self.cursor.fetchone()
            self.ra = result['ra0']
            self.dec = result['dec0']
            self.type = result['classification'] #TODO - update this to look at the type table and make a human readable type
            self.z = result['redshift']
            #These are set by the get_photometry method
            self.quality = {}
            self.jd = {}
            self.apparent_mag = {}
            self.apparent_mag_err = {}
            self.telescope = {}

    def get_sn_id(self):
        num_results = self.cursor.execute("SELECT DISTINCT targetid FROM targetnames WHERE name LIKE '%{}%'".format(self.name))
        result = self.cursor.fetchall()
        if num_results == 0:
            return None
        elif num_results>1:
            warnings.warn('Multiple targetids match your query. Please enter the targetid you would like to select')
            for iresult in result:
                target = self.cursor.execute("SELECT name FROM targetnames WHERE targetid = {}".format(iresult['targetid']))
                target_results = target.fetchall()
                print('targetid: {}'.format(iresult['targetid']))
                for itarget_result in target_results:
                    print(itarget_result['name'])
            targetid = input('Enter target id: ')
        else:
            targetid = result[0]['targetid']
        return targetid
 
    def get_photometry(self, band='all'):
        '''
        Get apparent magnitude of a specific band
        '''
        filter_dict =  define_filters()
        if band == 'all':
            self.cursor.execute("SELECT DISTINCT BINARY(filter) FROM photlco WHERE targetid={}".format(self.id))
            results = self.cursor.fetchall()
            bands = [iband['BINARY(filter)'].decode('utf-8') for iband in results]

        elif len(band)==1 or isinstance(band, str): #Single filter
            bands = [band]
        else: #multiple filters
            bands = band
        for ifilter in bands:   
            #self.cursor.execute("SELECT jd, mag, magerr, datatype, configuration  FROM photometry JOIN photometrysource on photometry.source=photometrysource.source WHERE targetid={} AND filter=BINARY('{}')".format(self.id, ifilter))
            sql_str = "SELECT mjd, mag, dmag, telescope, quality FROM photlco WHERE"\
                                + " wcs != 9999 AND (psf != 'X' AND psf != 9999) AND (zcat != 'X' OR abscat !='X') AND mag != 9999"\
                                + " AND targetid={} AND filter=BINARY('{}')".format(self.id, ifilter)
            self.cursor.execute(sql_str)
            results = self.cursor.fetchall()
            jd = []
            mag = []
            mag_err = []
            quality = []
            telescope = []
            for irow in results:
                if irow['dmag'] is None:
                    irow['dmag'] = 0
                jd.append(Time(float(irow['mjd']), format='mjd').jd)
                mag.append(float(irow['mag']))
                mag_err.append(float(irow['dmag']))
                quality.append(int(irow['quality']))
                telescope.append(irow['telescope'])
            self.jd[ifilter] = np.array(jd)
            self.apparent_mag[ifilter] = np.array(mag)
            self.apparent_mag_err[ifilter] = np.array(mag_err)
            self.quality[ifilter] = np.array(quality)
            self.telescope[ifilter] = np.array(telescope)
################

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
            the difference in days between the observation and the closest photometric point
    '''
    date_indx = np.argmin(np.abs(jd - date_obs.jd))
    closest_mag = phot_mag[date_indx]
    date_sep = (jd - date_obs.jd)[date_indx]
    return closest_mag, date_sep

def get_interpolated_photometry(date_obs, jd, phot_mag):
    '''
    Interpolate between photometric points in a given band to get the value at a single date
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
        interp_mag: float
            the magnitude interpolated to the date_obs given
        date_sep: tup
            the difference in days between the observation and the two points used in the interpolation
    '''

    before_indx = jd < date_obs.jd
    after_indx = jd >= date_obs.jd
    if (before_indx==False).all():
        before_sep = None
    else:
        before_sep = Time(jd[before_indx][-1], format='jd') - date_obs
    if (after_indx == False).all():
        after_sep = None
    else:
        after_sep =  Time(jd[after_indx][0],   format='jd') - date_obs
    if (before_sep is not None) and (after_sep is not None):
        interp_mag = np.interp(date_obs.jd, jd, phot_mag)
    else:
        interp_mag = None
    return interp_mag, (before_sep, after_sep)
    
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
                        header_date=False, sndavis=True, lightcurve=None, verbose=False):
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
    if lightcurve is None:
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
        if ifilter not in ['us', 'vs', 'bs', 'uw1', 'uw2', 'um1', 'um2', 'UVW1', 'UVW2', 'UVM1', 'UVM2']: #avoid swift filters
            cenwave = get_cenwave(ifilter)
            if (cenwave > 3500):
                mag, date_sep = get_interpolated_photometry(date_obs, lightcurve.jd[ifilter], lightcurve.apparent_mag[ifilter])
                if mag is not None:
                    print('using {}={} interpolated to {} from {} and {}'.format(ifilter, mag, date_obs.iso,
                                                                               (date_sep[0]+date_obs).iso,
                                                                               (date_sep[1]+date_obs).iso))
                    if ifilter.endswith('p'): #strip off p for LCO filters since qubascalespectra assumes filters are single character
                        band = band+ifilter[0]
                    else:
                        band=band+ifilter
                    mphot.append(mag)
                else:
                    if date_sep[0] is None:
                        print('No data for {} before {}'.format(ifilter, date_obs))
                    if date_sep[1] is None:
                        print('No data for {} after {}'.format(ifilter, date_obs))
    if len(mphot) > 0:
        if verbose:
            print('bands: {}, photometry: {}'.format(band, mphot))
        qubascalespectra.scale_spectrum(filename, band, mphot, filter_dir)
        shutil.move('log.txt', filename.split('.')[0]+'.log')
        with open(filename.split('.')[0]+'.log', 'a') as ofile:
            ofile.write(band+'\n')
            mphot_str = [str(iphot) for iphot in mphot]
            ofile.write(','.join(mphot_str)+'\n')
    else:
        print('WARNING: no photometry within {} days of observation date ({}) ({})'.format(max_sep, date_obs, filename))
    
        
            
    
        
