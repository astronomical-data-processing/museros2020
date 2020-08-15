# %matplotlib inline

import os
import sys
import numpy
import matplotlib
from rascil.processing_components import create_named_configuration
import argparse

# from matplotlib import plt.savefig
from astropy.coordinates import EarthLocation, SkyCoord, ITRS

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
from astropy.coordinates import get_body_barycentric, get_body, get_moon

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
from muser.data_models.muser_data import MuserData
from muser.data_models.muser_phase import MuserPhase
from muser.data_models.parameters import muser_path, muser_data_path
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


# def read_muser_data(file_name: str=''):


def main(args):
    if args.muser == '1':
        muser_array = 'MUSER1'
    else:
        muser_array = 'MUSER2'
    start_time = args.start
    end_time = args.end

    location = EarthLocation(lon=115.2505 * u.deg, lat=42.211833333 * u.deg, height=1365.0 * u.m)

    muser = MuserData(sub_array=1)
    if not muser.init_data_environment():
        print("No data environment prepared, exit.")
        return -1
    if not muser.search_first_file(frame_time=args.start):
        print("Cannot find observational data or not a MUSER file.")
        return -1
    if not muser.search_frame(args.start):
        print("Cannot locate the specified frame")
        return -1

    data_file_name = muser.current_file_name
    print("Checking MUSER File Information V20200801")
    print("First Observational Time {}".format(muser.current_frame_time.isot))
    # Check data
    print("Filename {} is a valid MUSER Data File.".format(data_file_name))
    print("Current Observational Time {}".format(muser.current_frame_utc_time.isot))
    print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser.is_loop_mode else "Non Loop", muser.frequency))
    print("Sub Band: {} - Sub Channel {}".format(muser.sub_band, muser.sub_channels))

    # count total frames
    total_frames = muser.count_frame_number(args.start, args.end)


    if muser.is_loop_mode:
        block_full_data = numpy.zeros([1, muser.antennas, muser.antennas, muser.sub_channels * muser.frame_number, 2],
                                      dtype='complex')
    else:
        block_full_data = numpy.zeros([1, muser.antennas, muser.antennas, muser.sub_channels, 1], dtype='complex')

    # Load Phase Calibration Data
    print("Loading Phase Calibration File")
    phase_cal = MuserPhase(muser.sub_array, muser.is_loop_mode, muser.current_frame_time)
    if not phase_cal.load_calibration_data():
        print("Cannot find phase calibration file. ")
        exit(1)
    print("File shape", phase_cal.phase_data.shape)

    # Load Phase Calibration Data Generated by 'export_phase_cal.py'

    # Create configuration of RASCIL
    muser_core = create_configuration(muser_array)

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
            for i in range(33 * 16):
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
        if not muser.read_full_frame():
            print("File reading error. ")
            exit(1)
        obs_time = muser.first_frame_time + 0.025 * u.second
        # TODO - J2000.0
        phasecentre = get_body('sun', obs_time, location)

        # convert phasecentre into ITRS coordinate
        c_ITRS = phasecentre.transform_to(ITRS(obstime=muser.current_frame_utc_time))
        local_ha = location.lon - c_ITRS.spherical.lon
        local_ha.wrap_at(24 * u.hourangle, inplace=True)

        print("Local Hour Angle: {}".format(local_ha.to('rad').value))
        times.append([local_ha.to('rad').value])
        integration_time.append(0.025)

        # phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')

        # visshape = [ntimes, nants, nants, nchan, npol]
        utc_time = Time('%04d-%02d-%2dT00:00:00' % (
            muser.current_frame_utc_time.datetime.year, muser.current_frame_utc_time.datetime.month,
            muser.current_frame_utc_time.datetime.day), format='isot')
        # Phase Calibration
        muser.phase_calibration(
            phase_cal.phase_data[..., muser.sub_band * muser.sub_channels:muser.sub_band * muser.sub_channels + 16,
            muser.polarization])

        # Inject data into blockvisibility
        vis_data[count, :, :, muser.sub_band * muser.sub_channels:muser.sub_band * muser.sub_channels + 16,
        1 - muser.polarization] = muser.block_data[:, :, :]
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
    bvis.vis[...] = vis_data[...]
    vis_list = []
    vis_list.append(bvis)

    export_file_name = muser_data_path(data_file_name) + '.ms'
    export_blockvisibility_to_ms(export_file_name, vis_list, source_name='SUN')

    # matplotlib.use('Agg')

    from matplotlib import pylab

    from matplotlib import pyplot as plt
    vt = convert_blockvisibility_to_visibility(bvis)
    plt.clf()
    plt.plot(vt.data['uvw'][:, 0], vt.data['uvw'][:, 1], '.', color='b')
    plt.plot(-vt.data['uvw'][:, 0], -vt.data['uvw'][:, 1], '.', color='r')
    plt.xlabel('U (wavelengths)')
    plt.ylabel('V (wavelengths)')
    plt.title("UV coverage")
    # plt.savefig(storedir + '/UV_coverage.pdf', format='pdf')
    plt.show()

    #
    # advice = advise_wide_field(vt, guard_band_image=3.0, delA=0.1, facets=1, wprojection_planes=1,
    #                            oversampling_synthesised_beam=4.0)
    # cellsize = advice['cellsize']

    # uvdist = numpy.sqrt(vt.data['uvw'][:, 0] ** 2 + vt.data['uvw'][:, 1] ** 2)
    #
    # model = create_image_from_visibility(vt, cellsize=cellsize, npixel=512)
    # dirty, sumwt = invert_list_serial_workflow([vt], [model], context='2d')[0]
    # psf, sumwt = invert_list_serial_workflow([vt], [model], context='2d', dopsf=True)[0]
    #
    # show_image(dirty)
    # print("Max, min in dirty image = %.6f, %.6f, sumwt = %f" % (dirty.data.max(), dirty.data.min(), sumwt))
    #
    # print("Max, min in PSF         = %.6f, %.6f, sumwt = %f" % (psf.data.max(), psf.data.min(), sumwt))
    # results_dir="/Users/f.wang"
    # export_image_to_fits(dirty, '%s/imaging_dirty.fits' % (results_dir))
    # export_image_to_fits(psf, '%s/imaging_psf.fits' % (results_dir))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='List Muser Data Information for Each Frame')
    parser.add_argument('-m', "--muser", type=str, default='1', help='The MUSER array')
    parser.add_argument('-c', "--calib", type=str, default='', help='The Calibration file name')
    parser.add_argument('-f', "--file", type=str, default='', help='The file name')
    parser.add_argument('-l', "--line", type=int, default=1, help='The number of frames')
    parser.add_argument('-s', "--start", type=str, default='', help='The beginning time ')
    parser.add_argument('-e', "--end", type=str, default='', help='The end time ')

    main(parser.parse_args())
    print('Done')
