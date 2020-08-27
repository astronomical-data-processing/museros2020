""" Coordinate scripts

"""

import numpy
import math
from astropy.units import u

def ecef_to_enu( x,  y,  z, lat0, lon0, h0):
    #  Convert to radians in notation consistent with the paper:
    a = 6378137.0         # WGS-84 Earth semimajor axis (m)

    b = 6356752.314245     # Derived Earth semiminor axis (m)
    f = (a - b) / a           # Ellipsoid Flatness
    f_inv = 1.0 / f       #Inverse flattening

    a_sq = a * a
    b_sq = b * b
    e_sq = f * (2 - f)    # Square of Eccentricity

    n_lambda = lat0*numpy.pi/180.
    phi = lon0*numpy.pi/180.
    s = numpy.sin(n_lambda)
    N = a / numpy.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(n_lambda)
    cos_lambda = math.cos(n_lambda)
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    xd = x - x0
    yd = y - y0
    zd = z - z0

    # This is the matrix multiplication
    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd
    return xEast,yNorth,zUp


def locxyz2itrf(lat, longitude, locx=0.0, locy=0.0, locz=0.0):
    """
    Returns the nominal ITRF (X, Y, Z) coordinates (m) for a point at "local"
    (x, y, z) (m) measured at geodetic latitude lat and longitude longitude
    (degrees).  The ITRF frame used is not the official ITRF, just a right
    handed Cartesian system with X going through 0 latitude and 0 longitude,
    and Z going through the north pole.  The "local" (x, y, z) are measured
    relative to the closest point to (lat, longitude) on the WGS84 reference
    ellipsoid, with z normal to the ellipsoid and y pointing north.
    """
    # from Rob Reid;  need to generalize to use any datum...
    import math
    phi, lmbda = map(math.radians, (lat, longitude))
    sphi = math.sin(phi)
    a = 6378137.0  # WGS84 equatorial semimajor axis
    b = 6356752.3142  # WGS84 polar semimajor axis
    ae = math.acos(b / a)
    N = a / math.sqrt(1.0 - (math.sin(ae) * sphi) ** 2)

    # Now you see the connection between the Old Ones and Antarctica...
    # Nploczcphimlocysphi = (N + locz) * pl.cos(phi) - locy * sphi
    Nploczcphimlocysphi = (N + locz) * math.cos(phi) - locy * sphi

    clmb = numpy.cos(lmbda)
    slmb = numpy.sin(lmbda)

    x = Nploczcphimlocysphi * clmb - locx * slmb
    y = Nploczcphimlocysphi * slmb + locx * clmb
    z = (N * (b / a) ** 2 + locz) * sphi + locy * math.cos(phi)

    return x, y, z


def compute_uvw(xyz, H, d, lat):
    """ Converts X-Y-Z coordinates into U-V-W
    Uses the transform from Thompson Moran Swenson (4.1, pg86)
    Parameters
    ----------
    xyz: should be a numpy array [x,y,z]
    H: float (degrees)
      is the hour angle of the phase reference position
    d: float (degrees)
      is the declination
    """
    sin = numpy.sin
    cos = numpy.cos
    trans1 = numpy.array([[0, -numpy.sin(lat), numpy.cos(lat)],
                    [1,  0,               0],
                    [0,  numpy.cos(lat), numpy.sin(lat)]])
    trans2 = numpy.array([
        [sin(H), cos(H), 0],
        [-sin(d) * cos(H), sin(d) * sin(H), cos(d)],
        [cos(d) * cos(H), -cos(d) * sin(H), sin(d)]
    ])
    uvw1 = xyz @ trans1.T
    uvw = uvw1 @ trans2.T
    return uvw