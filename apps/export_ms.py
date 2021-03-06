#! /usr/bin/env python

# %matplotlib inline

import os
import sys
import numpy
import copy
import matplotlib
from rascil.processing_components import create_named_configuration
import argparse
import logging


# from matplotlib import plt.savefig
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon

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

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

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


def main(args):
    log = logging.getLogger('muser')
    if len(args.log)>0:
        log.setLevel(level=logging.DEBUG)
        logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                            level=logging.DEBUG)
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)

        fh = logging.FileHandler(args.log,mode='w')
        fh.setLevel(logging.DEBUG)

        # 为logger对象添加句柄
        log.addHandler(console)
        log.addHandler(fh)

    if args.muser == 1:
        muser_array = 'MUSER1'
    else:
        muser_array = 'MUSER2'
    start_time = args.start
    end_time = args.end
    fringe = args.fringe
    if args.nolimit==0:
        nolimit = None
    else:
        nolimit = args.nolimit

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
    muser.search_frame(search_time=start_time)
    total_frames = muser.count_frame_number(start_time, end_time)
    print("Total {} frames will be processed.".format(total_frames))

    # Load Phase Calibration Data
    print("Loading Phase Calibration File")
    phase_cal = MuserPhase(muser.sub_array, muser.is_loop_mode, muser.current_frame_time, args.calib)
    if not phase_cal.load_calibration_data(file_name=args.calib):
        print("Cannot find phase calibration file. ")
        exit(1)

    # # Create configuration of RASCIL
    # xx,yy,zz = locxyz2itrf(42.211833333,115.2505,0,0,1365)
    muser_core = create_configuration(muser_array)
    # for x,y,z in muser_core.xyz:
    #     x,y,z = locxyz2itrf(42.211833333,115.2505,x,y,z+1365)
    #     print('{},{},{}'.format(x-xx,y-yy,z-zz))

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
    utc_times = []
    # Re-Search file
    if not muser.search_first_file(frame_time=args.start):
        print("Cannot find observational data or not a MUSER file.")
        return -1
    if not muser.search_frame(search_time=start_time):
        print("Cannot locate the specified frame")
        return -1
    log.info("Search file : {}".format(start_time))
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
            muser.delay_process('sun')

        if muser.is_loop_mode:
            obs_time = muser.first_frame_utc_time + 0.0125 * u.second
        else:
            obs_time = muser.first_frame_utc_time + 0.0015625 * u.second

        utc_times.append(obs_time)
        print("No.{} : Observation time (UTC) {}".format(count, obs_time))
        # Compute the position of the Sun
        Alpha, Delta, ha, Thete_z, Phi = get_sun(obs_time)

        times.append([(ha*u.deg).to('rad').value]) #[local_ha.to('rad').value])
        integration_time.append(0.025)

        phasecentre = SkyCoord(ra=Alpha * u.deg, dec=Delta * u.deg, frame='icrs', equinox='J2000')
        # visshape = [ntimes, nants, nants, nchan, npol]
        utc_time = Time('%04d-%02d-%02dT00:00:00' % (
            muser.current_frame_utc_time.datetime.year, muser.current_frame_utc_time.datetime.month,
            muser.current_frame_utc_time.datetime.day), format='isot')
        # Phase Calibration

        if args.calib == '':
            muser.phase_calibration(cal=phase_cal.phase_data)
        else:
            muser.phase_calibration(phai_sat=phase_cal.phase_data)

        # Inject data into blockvisibility
        if muser.is_loop_mode:
            vis_data[count, :, :, :, 1] = deepcopy(muser.block_full_data[:, :, :, 0])
            vis_data[count, :, :, :, 0] = deepcopy(muser.block_full_data[:, :, :, 1])
        else:
            vis_data[count, :, :, :, 0] = deepcopy(muser.block_full_data[:, :, :,0])

        count = count + 1

    times = numpy.array(times)
    integration_time = numpy.array(integration_time)
    if muser.is_loop_mode:
        bvis = create_blockvisibility(muser_core, times, frequency, phasecentre=phasecentre,
                                      weight=1.0, polarisation_frame=PolarisationFrame('circularnp'),
                                      channel_bandwidth=channelbandwidth,
                                      integration_time=integration_time,
                                      source='SUN',
                                      elevation_limit=nolimit,
                                      utc_time=utc_times)
    else:
        bvis = create_blockvisibility(muser_core, times, frequency, phasecentre=phasecentre,
                                      weight=1.0, polarisation_frame=PolarisationFrame('stokesI'),
                                      channel_bandwidth=channelbandwidth,
                                      integration_time=integration_time,
                                      source='SUN',
                                      elevation_limit=nolimit,
                                      utc_time=utc_times)
    bvis.data['vis'] = copy.deepcopy(vis_data)
    bvis.vis[...] = copy.deepcopy(vis_data[...])
    vis_list = []
    vis_list.append(bvis)

    # Output results
    if len(args.output)==0:
        output_time = Time(start_time, format='isot').datetime
        file_name = 'CSRH_%04d%02d%02d-%02d%02d%02d'%(output_time.year, output_time.month, output_time.day, output_time.hour, output_time.minute, output_time.second)
        export_file_name = muser_output_path(file_name) + '.ms'   #data_file_name
    else:
        export_file_name = muser_output_path(args.output) + '.ms'
    export_blockvisibility_to_ms(export_file_name, vis_list, source_name='SUN')
    print("Export file: {}".format(export_file_name))
    print("Done. ")

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Output Measurement Set Files')
    parser.add_argument('-m', "--muser", type=int, default=1, help='The MUSER array')
    parser.add_argument('-c', "--calib", type=str, default='', help='The Calibration file name')
    parser.add_argument('-f', "--file", type=str, default='', help='The file name')
    parser.add_argument('-s', "--start", type=str, default='', help='The beginning time ')
    parser.add_argument('-e', "--end", type=str, default='', help='The end time ')
    parser.add_argument('-t', "--fringe", type=str2bool, nargs = '?', const=True, default=False, help='Fringe Stop')
    parser.add_argument('-o', "--output", type=str, default='', help='The output file name')
    parser.add_argument('-l', "--log", type=str, default='', help='The output log file name')
    parser.add_argument('-n', "--nolimit", type=int, default=15, help='No limitation for the elevation')

    main(parser.parse_args())
