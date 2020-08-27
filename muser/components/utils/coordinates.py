""" Coordinate scripts

"""

import numpy
import math

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