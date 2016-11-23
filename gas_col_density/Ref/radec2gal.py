# from astropy.coordinates import ICRS, Galactic
# from astropy             import units as u

# c = ICRS(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree))

# #c = ICRS('00h42m44.3s +41d16m9s')
# print c.galactic

from astropy import units as u
from astropy.coordinates import SkyCoord

#c_icrs = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, frame='icrs')

c_icrs = SkyCoord('04h07m25.80s', '+03d40m10.4s', frame='icrs')
print c_icrs.galactic

c_icrs = SkyCoord('22h53m34.70s', '+16d11m25.70s', frame='icrs')
print c_icrs.galactic