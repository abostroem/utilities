import os

#Anaconda
import numpy as np
from astropy.time import Time, TimeDelta
from astropy import units as U
from astropy.io import ascii as asc

#Pip
import extinction
import pymysql
#Local
from define_filters import define_filters

CONFIG_DIR = os.environ['SNDAVIS']

class LightCurve(object):
    '''
    Caveats:
    * Assume plateau length is 100 days unless otherwise specified
    * Assume change in magnitude from peak to bottom of drop at end of plateau is 5.5 mag
    * I'm calculating the age from the peak not from explosion
    '''
    def __init__(self, peak_date, peak_mag, plateau_length = 100, plateau_delta_mag = 5.5):
        self.peak_time = Time(peak_date, format='iso')
        self.peak_mag = peak_mag
        self.plateau_length = TimeDelta(plateau_length*U.day)
        self.plateau_delta_mag = plateau_delta_mag
    def calc_mag_at_first_obs(self, first_obs_date):
        self.drop_date = self.peak_time+self.plateau_length
        self.cobalt_decay_delta_t = Time(first_obs_date, format='iso') - self.drop_date
        self.cobalt_decay_delta_mag = self.cobalt_decay_delta_t/(100.*U.day)
        self.mag_at_first_obs = self.peak_mag + self.plateau_delta_mag + self.cobalt_decay_delta_mag
        self.age_at_first_obs = self.plateau_length + self.cobalt_decay_delta_t
        print('Mag on {} = {}'.format(first_obs_date, self.mag_at_first_obs))
        print('Age on {} = {}'.format(first_obs_date, self.age_at_first_obs))
        
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
    abs_mag_err = np.sqrt(app_mag_err**2 + dist_mod_err**2)
    return abs_mag, abs_mag_err

def database_config(config_file):
    config = asc.read(config_file, data_start=0)
    config_dict = {}
    for keyword, arg in config:
        config_dict[keyword]=arg
    return config_dict

class LightCurve2(object):
    '''
    Retrieve information from sndavis database
    '''
    def __init__(self, name):
        self.name=name
        config_dict = database_config(os.path.join(CONFIG_DIR, 'configure'))
        if 'port' in config_dict.keys():
            self.db = pymysql.connect(user=config_dict['mysqluser'], 
                         host=config_dict['hostname'], 
                         password=config_dict['mysqlpasswd'],
                         database=config_dict['database'],
                        port=int(config_dict['port']))
        else:
            self.db = pymysql.connect(user=config_dict['mysqluser'], 
                         host=config_dict['hostname'], 
                         password=config_dict['mysqlpasswd'],
                         database=config_dict['database'])
        self.cursor = self.db.cursor()
        self.id = self.get_sn_id()
        self.cursor.execute('SELECT * FROM idsupernovae WHERE id = {}'.format(self.id))
        result = self.cursor.fetchone()
        self.ra = result[1]
        self.dec = result[2]
        self.type = result[3] #TODO - update this to look at the type table and make a human readable type
        self.dist_mod = result[10]
        self.dist_mod_err = result[11]
        self.ebv_mw = result[13]
        self.ebv_mw_err = result[14]
        self.ebv_host = result[15] 
        self.ebv_host_err = result[16]
        self.quality = result[17]
        self.jdexpl = result[18]
        self.jdexpl_err = result[19]
        self.jd = {}
        self.apparent_mag = {}
        self.apparent_mag_err = {}
        self.phase = {}
        self.abs_mag = {}
        self.abs_mag_err = {}
        self.slopes = {'s1':{}, 's2':{}, 's50':{}, 'tail':{}, 
                        's1_range':{}, 's2_range':{}, 's50_range':{}, 'tail_range':{},
                        's1_err':{}, 's2_err':{}, 's50_err':{}, 'tail_err':{}}
        #TODO add an absolute magnitude error
             
    def get_sn_id(self):
        self.cursor.execute('SELECT targetid FROM supernovanames WHERE name = {}'.format('"{}"'.format(self.name)))
        result = self.cursor.fetchone()
        return result[0]
 
    def get_photometry(self, band='all'):
        '''
        Get apparent magnitude of a specific band
        '''
        if band == 'all':
            self.cursor.execute("SELECT DISTINCT BINARY(filter) FROM photometry WHERE targetid={}".format(self.id))
            results = self.cursor.fetchall()
            bands = [str(iband[0]).strip('b').strip("'") for iband in results]

        elif len(band)==1: #Single filter
            bands = [band]
        else: #multiple filters
            bands = band
        for ifilter in bands:   
            self.cursor.execute("SELECT jd, mag, magerr  FROM photometry WHERE targetid={} AND BINARY filter='{}'".format(self.id, ifilter))
            results = self.cursor.fetchall()
            jd = []
            mag = []
            mag_err = []
            for irow in results:
                jd.append(float(irow[0]))
                mag.append(float(irow[1]))
                mag_err.append(float(irow[2]))
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
                                                         BINARY filter="{}"'.format(self.id, slope_type, band))
        result = self.cursor.fetchone()
        if slope_type == 's1':
            self.slopes['s1'][band] = result[2]
            self.slopes['s1_err'][band] = result[3]
            self.slopes['s1_range'][band] = (result[4], result[5])
        elif slope_type == 's2':
            self.slopes['s2'][band] = result[2]
            self.slopes['s2_err'][band] = result[3]
            self.slopes['s2_range'][band] = (result[4], result[5])
        elif slope_type == 's50':
            self.slopes['s50'][band] = result[2]
            self.slopes['s50_err'][band]= result[3]
            self.slopes['s50_range'][band] = (result[4], result[5])
        elif slope_type == 'tail':
            self.slopes['tail'][band] = result[2]
            self.slopes['tail_err'][band] = result[3]
            self.slopes['tail_range'][band] = (result[4], result[5])
        
    
