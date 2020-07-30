"""Various constants required by Skyfield."""

from numpy import array

# Angles.
ASEC360 = 1296000.0
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321
TAU = 6.283185307179586476925287
tau = TAU  # lower case, for symmetry with math.pi

# Physics.
C_AUDAY = 173.1446326846693
C = 299792458.0

# Earth and its orbit.
ANGVEL = 7.2921150e-5
AU_M = 1.4959787069098932e+11
AU_KM = 1.4959787069098932e+8
ERAD = 6378136.6
IERS_2010_INVERSE_EARTH_FLATTENING = 298.25642

PSI_COR = 0.0
EPS_COR = 0.0

# Heliocentric gravitational constant in meters^3 / second^2, from DE-405.
GS = 1.32712440017987e+20

# Time.
T0 = 2451545.0
B1950 = 2433282.4235
DAY_S = 86400.0

# import numpy
# numpy.set_printoptions(formatter={'float': repr})
# from .constants import ASEC2RAD, T0
# from .nutationlib import mean_obliquity
# from .functions import rot_x
# ecliptic_obliquity_radians = mean_obliquity(T0) * ASEC2RAD
# print(repr(rot_x(-ecliptic_obliquity_radians)))

rotation_to_ecliptic = array(((1.0, 0.0, 0.0),
       (0.0, 0.91748214306524178, 0.39777696911260602),
       (0.0, -0.39777696911260602, 0.91748214306524178)))
