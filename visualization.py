from astropy.visualization import ZScaleInterval
def zscale(img):
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(img)
    return vmin, vmax