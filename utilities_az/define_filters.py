import collections
def define_filters():
    # Filter parameters
    # Zero point vegamag computed with synphot
    bandpar=collections.OrderedDict()
    band_label= ['bandw','fwhm','avgwv','equvw','zpoint']
    bandpar['A']=[348.43,820.5,2650.6,770.45,3.818e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['uw1']=[348.43,820.5,2650.6,770.45,3.818e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2
    bandpar['UVW1']=[348.43,820.5,2650.6,770.45,3.818e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['S']=[285.12,671.4,2136.7,639.45,4.825e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['uw2']=[285.12,671.4,2136.7,639.45,4.825e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['UVW2']=[285.12,671.4,2136.7,639.45,4.825e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['D']=[189.46,446.14,2269.2,519.03,4.321e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['um2']=[189.46,446.14,2269.2,519.03,4.321e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2
    bandpar['UVM2']=[189.46,446.14,2269.2,519.03,4.321e-9] # SWIFT A-> UW1, S -> UW2, D -> UM2 
    bandpar['us'] = [233.63, 550.15, 3442.6, 37007, 3.605e-9] # swift filter through andpar in synphot
    bandpar['u']=[194.41,457.79, 3561.8, 60.587,3.622e-9] # sloan
    bandpar['U']=[205.79,484.6,3652,542.62,4.327e-9] #landolt & johnson buser and kurucz 78
    bandpar['bs'] = [291.93, 687.45, 4438.7, 32718, 6.359e-9] # swift filter through andpar in synphot
    bandpar['B']=[352.94,831.11,4448,1010.3,6.09e-9] #landolt & johnson buser and kurucz 78
    bandpar['g']=[394.17,928.19, 4718.9, 418.52,5.414e-9] # sloan
    bandpar['gp']=[394.17,928.19, 4718.9, 418.52,5.414e-9] # LCO g band
    bandpar['vs'] = [232.12, 546.61, 5436, 15447, 3.657e-9] # swift filter through andpar in synphot
    bandpar['V']=[351.01,826.57,5505,870.65,3.53e-9] #landolt & johnson buser and kurucz 78
    bandpar['r']=[344.67,811.65, 6185.2, 546.14,2.472e-9] # sloan
    bandpar['rp']=[344.67,811.65, 6185.2, 546.14,2.472e-9] # LCO r band
    bandpar['Open']=[344.67,811.65, 6185.2, 546.14,2.472e-9] # DLT40 Prompt5
    bandpar['Clear']=[344.67,811.65, 6185.2, 546.14,2.472e-9] # DLT40 Meckering
    bandpar['R']=[589.71,1388.7,6555,1452.2,2.104e-9] #landolt & cousin bessel 83
    bandpar['i']=[379.57,893.82, 7499.8, 442.15,1.383e-9] # sloan
    bandpar['ip']=[379.57,893.82, 7499.8, 442.15,1.383e-9] # LCO i band
    bandpar['I']=[381.67,898.77,7900.4,1226,1.157e-9] #landolt & cousin bessel 83
    bandpar['z']=[502.45,1183.2, 8961.5, 88.505,8.15e-10] # sloan
    bandpar['Y']=[747.1,1759.3,10370,2034,3.05e-10] #bessel in bessel and brett 88
    bandpar['J']=[747.1,1759.3,12370,2034,3.05e-10] #bessel in bessel and brett 88
    bandpar['H']=[866.55,2040.6,16471,2882.8,1.11e-10] #bessel in bessel and brett 88
    bandpar['K']=[1188.9,2799.6,22126,3664.3,3.83e-11] #bessel in bessel and brett 88
    bandpar['L']=[0.,9000.,34000,0.,8.1e-12] # ASIAGO PHOTOMETRIC DATABASE
    bandpar['M']=[0.,11000.,50000,0.,2.2e-12] # ASIAGO PHOTOMETRIC DATABASE
    bandpar['N']=[0.,60000.,102000,0.,1.23e-13] # ASIAGO PHOTOMETRIC DATABASE
    bandpar['wfc3_f814w']=[672.53, 1583.68, 8081.89, 367.79] #Calculated with PySynphot ObsBandpass('wfc3,uvis1,f814w')
    return bandpar

def get_cenwave(band):
    band_dict = define_filters()
    if band in band_dict.keys():
        band_info = band_dict[band]
        cenwave = band_info[2]
        return cenwave
    else:
        print('band {} not in band dict keys: {}'.format(band, band_dict.keys()))
        raise 'KeyError' #TODO: figure out why this doesn't work