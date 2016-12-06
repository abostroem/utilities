from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib import pyplot
import matplotlib
import numpy as np


def zscale(img):
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(img)
    return vmin, vmax
    
def make_color_wheel(obj_list, cmap='rainbow'):
    cmap = matplotlib.cm.get_cmap('gist_rainbow')
    obj_color_list = cmap(np.linspace(0, 1, len(obj_list)))
    return obj_color_list
    
    
def plot_multiple_spectra(file_list, xcol='WAVELENGTH', ycol='FLUX', ax=None, colors = [], separation=0.5):
    if ax is None:
        fig = pyplot.figure(figsize=[10, 30])
        ax = fig.add_subplot(1,1,1)
    if len(colors)==0:
        colors = make_color_wheel(file_list)
    for indx, ifile in enumerate(file_list):
        ofile = fits.open(ifile)
        if indx ==0:
            scale1 = np.sum(ofile[1].data[ycol])
            offset = separation*np.max(ofile[1].data[ycol])
        scale2 = np.sum(ofile[1].data[ycol])
        l1=ax.plot(ofile[1].data[xcol], ofile[1].data[ycol]*scale1/scale2+offset*indx)
        if len(colors)>0:
            l1.set_color(colors[indx])
    return ax
    
