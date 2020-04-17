import os
from utilities_az import supernova

def test_LightCurveLCO():
    if os.path.exists('/Volumes/quince_external'):
        sn = supernova.LightCurveLCO('2018zd')                                                                                  
        sn.get_photometry()