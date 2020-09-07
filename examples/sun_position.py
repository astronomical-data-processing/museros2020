from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, ITRS
from astropy.coordinates import get_body_barycentric, get_body, get_moon, get_sun
import astropy.coordinates as coord

def sun_position(obs_time):
    location = EarthLocation(lon=115.2505 * u.deg, lat=42.211833333 * u.deg, height=1365.0 * u.m)

    solar_system_ephemeris.set('de432s')

    phasecentre = get_body('sun', obs_time, location, ephemeris='jpl')
    print('SUN: RA:{} Dec:{}'.format(phasecentre.ra.value, phasecentre.dec.value))
    # convert phasecentre into ITRS coordinate
    c_ITRS = phasecentre.transform_to(ITRS(obstime=obs_time))
    local_ha = location.lon - c_ITRS.spherical.lon
    local_ha.wrap_at(24 * u.hourangle, inplace=True)

    print("UTC: {} Local Hour Angle: {}".format(obs_time, local_ha.to('deg').value))
    obs_time1 = Time(obs_time,scale='utc', location=location)
    lst = obs_time1.sidereal_time('mean')
    local_ha1 = lst.to('deg').value - phasecentre.ra.value
    print("UTC: {} Local Hour Angle: {}".format(obs_time, local_ha1*u.deg))

    sun = coord.get_sun(obs_time)
    sun_icrs = sun.transform_to(coord.ICRS)
    print(sun_icrs)
def main(args):
    obs_time = Time(args.start, format='isot', scale='utc')
    sun_position(obs_time)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Output The SUN's position")
    parser.add_argument('-s', "--start", type=str, default='', help='The beginning time ')
    main(parser.parse_args())
