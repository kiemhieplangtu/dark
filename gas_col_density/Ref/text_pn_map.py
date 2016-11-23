import healpy as hp
import numpy as np
import pylab as pl

hp.mollview(title="Galactic map of Test Fields")
# hp.graticule()

x = [-37, 88, -137, -139, -136, -44]
y = [27, -60, -1.4, -50, -77, -46]
lab = ['TF0.1', 'TF0.2', 'TF0.3', 'TF0.4', 'TF0.5', 'TF0.6' ]

hp.projscatter(x, y, lonlat=True, coord='G')

hp.projtext(-37., 27., 'TF0.1', lonlat=True, coord='G')
hp.projtext(88, -60, 'TF0.2', lonlat=True, coord='G')
hp.projtext(-137, -1.4, 'TF0.3', lonlat=True, coord='G')
hp.projtext(-139, -50, 'TF0.4', lonlat=True, coord='G')
hp.projtext(-136, -77, 'TF0.5', lonlat=True, coord='G')
hp.projtext(-44, -46, 'TF0.6', lonlat=True, coord='G')

# equateur_lon = [-45.,45.]
# equateur_lat = [-30.,30.]
# hp.projplot(equateur_lon, equateur_lat, 'ro-', lonlat=True, coord='G')
pl.show()