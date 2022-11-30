import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from matplotlib.colors import PowerNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm

#WCS coordinate system
w = wcs.WCS(naxis=2)
w.wcs.crpix = [23.5, 23.5]
w.wcs.cdelt = np.array([-0.0035, 0.0035])
w.wcs.crval = [266.8451, -28.151658]
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w.wcs.set_pv([(2, 1, 45.0)])

#generate an array as image test
data = (np.arange(10000).reshape((100,100)))

#display image
fig = plt.figure()
ax = plt.gca(projection=w)
graf = ax.imshow(data, origin='lower', cmap=cm.viridis, norm=PowerNorm(1))

#colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%")
cbar = fig.colorbar(graf, cax=cax, orientation='horizontal')
cax.xaxis.set_ticks_position('top')
fig.show()
plt.show()