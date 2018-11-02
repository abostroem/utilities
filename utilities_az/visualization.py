from astropy.visualization import ZScaleInterval
from astropy.io import fits
from matplotlib import pyplot as plt
from cycler import cycler
import matplotlib as mpl
import numpy as np


def zscale(img):
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(img)
    return vmin, vmax
    
def make_color_wheel(obj_list, cmap='gist_rainbow'):
    cmap = mpl.cm.get_cmap(cmap)
    obj_color_list = cmap(np.linspace(0, 1, len(obj_list)))
    return obj_color_list
    
    
def plot_multiple_spectra(file_list, xcol='WAVELENGTH', ycol='FLUX', ax=None, colors = [], separation=0.5):
    if ax is None:
        fig = plt.figure(figsize=[10, 30])
        ax = fig.add_subplot(1,1,1)
    if len(colors)==0:
        colors = make_color_wheel(file_list)
    for indx, ifile in enumerate(file_list):
        ofile = fits.open(ifile)
        if indx ==0:
            scale1 = np.sum(ofile[1].data[ycol])
            offset = separation*np.max(ofile[1].data[ycol])
        scale2 = np.sum(ofile[1].data[ycol])
        l1=ax.plot(ofile[1].data[xcol], ofile[1].data[ycol]*scale1/scale2+offset*indx, color = colors[indx])
    return ax
    
def set_bright_colors():
    bright = ['#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377']
    plt.rc('axes', prop_cycle=(cycler('color', bright)))
    
def make_rainbow_cm():
    cols = [(0,0,0)]
    for x in np.linspace(0,1, 254):
        rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
        gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
        bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
        cols.append((rcol, gcol, bcol))

    cols.append((1,1,1))
    cm_rainbow = mpl.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
    return cm_rainbow