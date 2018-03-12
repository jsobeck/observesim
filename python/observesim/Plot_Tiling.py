import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import mpl_toolkits.basemap as basemap
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from astropy import units as u
from astropy.coordinates import Angle
%matplotlib 

##APO
apo     = Table.read('/Users/jsobeck/observesim/data/tiling/straw-apo.fits')
ra_apo  = coord.Angle(apo['racen']*u.degree)
ra_apo  = ra_apo.wrap_at(180*u.deg)
dec_apo = coord.Angle(apo['deccen']*u.degree)
fig     = plt.figure(figsize=(8,6))
ax      = fig.add_subplot(111, projection="aitoff")
ax.scatter(ra_apo.radian, dec_apo.radian, s=1, alpha=0.5, color = 'red')
plt.grid(True)
plt.savefig('SDSSV-APO_Tile_Centers.pdf')

##LCO
lco     = Table.read('/Users/jsobeck/observesim/data/tiling/straw-lco.fits')
ra_lco  = coord.Angle(lco['racen']*u.degree)
ra_lco  = ra_lco.wrap_at(180*u.deg)
dec_lco = coord.Angle(lco['deccen']*u.degree)
fig     = plt.figure(figsize=(8,6))
ax      = fig.add_subplot(111, projection="aitoff")
ax.scatter(ra_lco.radian, dec_lco.radian, s=1, alpha=0.5, color = 'navy')
plt.grid(True)
plt.savefig('SDSSV-LCO_Tile_Centers.pdf')

#################### M. Blanton #################
"""Plot ra and dec"""
m = basemap.Basemap(projection='moll', lon_0=270, resolution='c')

# draw parallels and meridians.
m.drawparallels(np.arange(-90., 120., 30.),
                    linewidth=0.5,
                    labels=[1, 0, 0, 0],
                    labelstyle='+/-')
m.drawmeridians(np.arange(0., 420., 60.), linewidth=0.5)
m.drawmapboundary()

radius = 1.5
boundary = 90.
ncirc = 60
theta = np.arange(ncirc) * np.pi * 2. / np.float32(ncirc - 1)

def convert_radec(m, ra, dec):
        return m(((360. - ra) + 180.) % 360., dec, inverse=False)

for (ra, dec) in zip(ra, dec):
    cra = (ra + radius / np.cos(dec * np.pi / 180.) * np.cos(theta))
    cdec = dec + radius * np.sin(theta)

    nbelow = np.count_nonzero(cra < boundary)
    if(nbelow > 0):
        ibelow = np.nonzero(cra < boundary)[0]
        (xx, yy) = convert_radec(m, cra[ibelow], cdec[ibelow])
        plt.plot(xx, yy, linewidth=0.25, color='blue')

    nabove = np.count_nonzero(cra > boundary)
    if(nabove > 0):
        iabove = np.nonzero(cra > boundary)[0]
        (xx, yy) = convert_radec(m, cra[iabove], cdec[iabove])
        plt.plot(xx, yy, linewidth=0.25, color='blue')

plt.show()
