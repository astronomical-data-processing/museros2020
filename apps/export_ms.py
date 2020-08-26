# %matplotlib inline

import os
import sys
import numpy
import copy
import matplotlib
from rascil.processing_components import create_named_configuration
import argparse

# from matplotlib import plt.savefig
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, AltAz

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# sys.path.append(os.path.join('..','..'))

from rascil.data_models.parameters import rascil_path

# results_dir = rascil_path('test_results')

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord

from rascil.data_models.polarisation import PolarisationFrame
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from muser.components.ephem.sun_position import get_sun
from copy import deepcopy
from rascil.processing_components import image_raster_iter
from rascil.processing_components import create_visibility
# from rascil.processing_components import sum_visibility
from rascil.processing_components import vis_timeslices, vis_wslices
from rascil.processing_components import create_configuration_from_file
from rascil.processing_components import create_skycomponent, find_skycomponents, \
    find_nearest_skycomponent, insert_skycomponent
from rascil.processing_components import show_image, export_image_to_fits, qa_image, smooth_image
from rascil.processing_components import advise_wide_field, create_image_from_visibility, \
    predict_skycomponent_visibility
from rascil.processing_components.util import azel_to_hadec
from muser.data_models.muser_data import MuserData
from muser.data_models.muser_phase import MuserPhase
from muser.data_models.parameters import muser_path, muser_data_path, muser_output_path

from rascil.processing_components.visibility.coalesce import convert_visibility_to_blockvisibility, \
    convert_blockvisibility_to_visibility
import logging

try:
    import casacore
    from casacore.tables import table  # pylint: disable=import-error
    from rascil.processing_components.visibility.base import create_blockvisibility, create_blockvisibility_from_ms
    from rascil.processing_components.visibility.base import export_blockvisibility_to_ms

    run_ms_tests = True
#            except ModuleNotFoundError:
except:
    run_ms_tests = False


def create_configuration(name: str = 'LOWBD2', **kwargs):
    from muser.data_models.parameters import muser_path
    import os
    conf_dir = muser_path('configurations')

    location = EarthLocation(lon=115.2505 * u.deg, lat=42.211833333 * u.deg, height=1365.0 * u.m)
    if name == 'MUSER2':
        antfile = os.path.join(conf_dir, 'muser-2.csv')
        lowcore = create_configuration_from_file(antfile=antfile,
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=2.0, name='MUSER', location=location, **kwargs)
    else:
        antfile = os.path.join(conf_dir, 'muser-1.csv')
        lowcore = create_configuration_from_file(antfile=antfile,
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=4.5, name='MUSER', location=location, **kwargs)
    return lowcore


def enu2xyz(lat, locx=0.0, locy=0.0, locz=0.0):
    pass

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
# def read_muser_data(file_name: str=''):


def main(args):
    if args.muser == 1:
        muser_array = 'MUSER1'
    else:
        muser_array = 'MUSER2'
    start_time = args.start
    end_time = args.end
    fringe = args.fringe

    location = EarthLocation(lon=115.2505 * u.deg, lat=42.211833333 * u.deg, height=1365.0 * u.m)

    muser = MuserData(sub_array=args.muser)
    if not muser.init_data_environment():
        print("No data environment prepared, exit.")
        return -1
    if not muser.search_first_file(frame_time=args.start):
        print("Cannot find observational data or not a MUSER file.")
        return -1
    data_file_name = os.path.basename(muser.current_file_name)
    print("Checking MUSER File Information V20200801")
    print("First Observational Time {}".format(muser.current_frame_time.isot))
    # Check data
    print("Filename {} is a valid MUSER Data File.".format(data_file_name))
    print("Current Observational Time {}".format(muser.current_frame_utc_time.isot))
    print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser.is_loop_mode else "Non Loop", muser.frequency))
    print("Sub Band: {} - Sub Channel {}".format(muser.sub_band, muser.sub_channels))

    # count total frames
    muser.search_frame(search_time=args.start)
    total_frames = muser.count_frame_number(args.start, args.end)

    # Load Phase Calibration Data
    print("Loading Phase Calibration File")
    phase_cal = MuserPhase(muser.sub_array, muser.is_loop_mode, muser.current_frame_time)
    if not phase_cal.load_calibration_data():
        print("Cannot find phase calibration file. ")
        exit(1)
    print("File shape", phase_cal.phase_data.shape)

    # # Create configuration of RASCIL
    # xx,yy,zz = locxyz2itrf(42.211833333,115.2505,0,0,1365)
    muser_core = create_configuration(muser_array)
    # for x,y,z in muser_core.xyz:
    #     x,y,z = locxyz2itrf(42.211833333,115.2505,x,y,z+1365)
    #     print('{},{},{}'.format(x-xx,y-yy,z-zz))

    # Create Phase Centre
    # TODO: need to compute the position of the SUN
    freq = []
    if muser.is_loop_mode:
        if muser.sub_array == 1:
            for i in range(64):
                freq.append(400e6 + 25e6 * i)
            channelbandwidth = numpy.array([25e6] * 64)
        else:
            for i in range(33 * 16):
                freq.append(2e9 + 25e6 * i)
            channelbandwidth = numpy.array([25e6] * 16 * 33)
    else:
        if muser.sub_array == 1:
            for i in range(16):
                freq.append(muser.frequency + 25e6 * i)
            channelbandwidth = numpy.array([25e6] * 16)
        else:
            for i in range(16):
                freq.append(muser.frequency + 25e6 * i)
            channelbandwidth = numpy.array([25e6] * 16)

    frequency = numpy.array(freq)
    integration_time = []  # numpy.array([0.025])
    times = []
    solar_system_ephemeris.set('de432s')

    # Re-Search file
    if not muser.search_first_file(frame_time=args.start):
        print("Cannot find observational data or not a MUSER file.")
        return -1
    if not muser.search_frame(args.start):
        print("Cannot locate the specified frame")
        return -1

    count = 0
    # total_frames = 1
    if muser.is_loop_mode:
        vis_data = numpy.zeros(
            (total_frames, muser.antennas, muser.antennas, muser.sub_channels * muser.frame_number, 2), dtype='complex')
    else:
        vis_data = numpy.zeros(
            (total_frames, muser.antennas, muser.antennas, 16, 1), dtype='complex')

    while count < total_frames:
        if not muser.read_full_frame(read_data=True):
            print("File reading error. ")
            exit(1)
        # Delay processing for the Sun
        if fringe:
            if muser.sub_array == 2:
                if muser.current_frame_header.strip_switch == 0xCCCCCCCC:
                    muser.delay_process("sun")
            else:
                muser.delay_process('sun')

        obs_time = muser.first_frame_utc_time #+ 0.0125 * u.second
        # TODO - J2000.0
        Alpha, Delta, ha, Thete_z, Phi = get_sun(obs_time)
        # with solar_system_ephemeris.set('de432s'):
        #     p1 = get_body('sun', muser.first_frame_utc_time, location)
        # phasecentre = get_body('sun', muser.first_frame_utc_time,ephemeris='de432s')
        # phasecentre = get_body('sun', muser.first_frame_utc_time,location=location,ephemeris='de432s')

        obs_time = Time(muser.first_frame_utc_time, location=location)
        lst_apparent = obs_time.sidereal_time('apparent')
        times.append([(ha*u.deg).to('rad').value]) #[local_ha.to('rad').value])
        integration_time.append(0.025)

        phasecentre = SkyCoord(ra=Alpha * u.deg, dec=Delta * u.deg, frame='icrs', equinox='J2000')

        # visshape = [ntimes, nants, nants, nchan, npol]
        utc_time = Time('%04d-%02d-%2dT00:00:00' % (
            muser.current_frame_utc_time.datetime.year, muser.current_frame_utc_time.datetime.month,
            muser.current_frame_utc_time.datetime.day), format='isot')
        # Phase Calibration
        muser.phase_calibration(phase_cal.phase_data)

        # Inject data into blockvisibility
        vis_data[count, :, :, :, 1] = deepcopy(muser.block_full_data[:, :, :, 0])
        vis_data[count, :, :, :, 0] = deepcopy(muser.block_full_data[:, :, :, 1])
        count = count + 1

    times = numpy.array(times)
    integration_time = numpy.array(integration_time)
    if muser.is_loop_mode:
        bvis = create_blockvisibility(muser_core, times, frequency, phasecentre=phasecentre,
                                      weight=1.0, polarisation_frame=PolarisationFrame('circularnp'),
                                      channel_bandwidth=channelbandwidth,
                                      integration_time=integration_time,
                                      source='SUN',
                                      utc_time=utc_time)
    else:
        bvis = create_blockvisibility(muser_core, times, frequency, phasecentre=phasecentre,
                                      weight=1.0, polarisation_frame=PolarisationFrame('stokesI'),
                                      channel_bandwidth=channelbandwidth,
                                      integration_time=integration_time,
                                      source='SUN',
                                      utc_time=utc_time)
    bvis.data['vis'] = copy.deepcopy(vis_data)
    bvis.vis[...] = copy.deepcopy(vis_data[...])
    vis_list = []
    vis_list.append(bvis)

    export_file_name = muser_output_path(data_file_name) + '.ms'
    export_blockvisibility_to_ms(export_file_name, vis_list, source_name='SUN')

    print("Done. ")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Output Measurement Set Files')
    parser.add_argument('-m', "--muser", type=int, default=1, help='The MUSER array')
    parser.add_argument('-c', "--calib", type=str, default='', help='The Calibration file name')
    parser.add_argument('-f', "--file", type=str, default='', help='The file name')
    parser.add_argument('-l', "--line", type=int, default=1, help='The number of frames')
    parser.add_argument('-s', "--start", type=str, default='', help='The beginning time ')
    parser.add_argument('-e', "--end", type=str, default='', help='The end time ')
    parser.add_argument('-t', "--fringe", type=bool, default=False, help='Fringe Stop')
    main(parser.parse_args())
