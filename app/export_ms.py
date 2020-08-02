# %matplotlib inline

import os
import sys
import numpy
import matplotlib
from rascil.processing_components import create_named_configuration
import argparse

# from matplotlib import plt.savefig
from astropy.coordinates import EarthLocation, SkyCoord

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# sys.path.append(os.path.join('..','..'))

from rascil.data_models.parameters import rascil_path
# results_dir = rascil_path('test_results')

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord

from rascil.data_models.polarisation import PolarisationFrame

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

    location = EarthLocation(lon=115.2505, lat=42.211833333, height=1365.0)
    if name == 'MUSER2':
        antfile = os.path.join(conf_dir, 'muser-2.csv')
        lowcore = create_configuration_from_file(antfile=antfile,
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=2.0, name='MUSER', location=location, **kwargs)
    else:
        antfile = os.path.join(conf_dir, 'muser-1.csv')
        lowcore = create_configuration_from_file(antfile=antfile,
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=4.0, name='MUSER', location=location, **kwargs)
    return lowcore

# def read_muser_data(file_name: str=''):


def main(args):

    if args.muser == '1':
        muser_array='MUSER1'
    else:
        muser_array='MUSER2'

    data_file_name = 'CSRH_20151122-093500_89058131'
    file_name = muser_data_path(data_file_name)
    muser = MuserData(sub_array=1, file_name=file_name)

    if not muser.check_muser_file():
        print("Cannot find observational data or not a MUSER file.")
        exit(1)
    print("Checking MUSER File Information V20200801")
    print("First Observational Time {}".format(muser.current_frame_time.isot))
    # Check data
    muser.search_frame('2015-11-22T09:41:00')
    print("Filename {} is a valid MUSER Data File.".format(file_name))
    print("Current Observational Time {}".format(muser.current_frame_time.isot))
    print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser.is_loop_mode else "Non Loop", muser.frequency))
    print("Sub Band: {} - Sub Channel {}".format(muser.sub_band,muser.sub_channels))
    muser.read_data()
    muser.delay_process('sun')
    # Create configuration
    muser_core = create_configuration(muser_array)

    # Create Phase Centre
    # TODO: need to compute the position of the SUN
    local_time = Time(muser.current_frame_utc_time, location=('115.2505d', '42.211833333d'))
    sidereal_time = local_time.sidereal_time('apparent')
    print("Sidereal time: {}".format(sidereal_time))
    times = numpy.array([sidereal_time.hour]) #muser.current_frame_time.datetime])
    freq = []
    for i in range(16):
        freq.append(muser.frequency+25e6*i)
    frequency = numpy.array(freq)
    channelbandwidth = numpy.array([25e6]*16)
    reffrequency = numpy.max(frequency)

    from astropy.coordinates import solar_system_ephemeris, EarthLocation
    from astropy.coordinates import get_body_barycentric, get_body, get_moon
    location = EarthLocation(lon=115.2505, lat=42.211833333, height=1365.0)
    solar_system_ephemeris.set('de432s')
    phasecentre = get_body('sun', muser.current_frame_utc_time, location)
    # phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')

    # visshape = [ntimes, nants, nants, nchan, npol]
    bvis = create_blockvisibility(muser_core, times, frequency, phasecentre=phasecentre,
                                  weight=1.0, polarisation_frame=PolarisationFrame('stokesI'),
                                  channel_bandwidth=channelbandwidth)
    bvis.vis[:,...,:] = muser.block_data
    vis_list = []
    vis_list.append(bvis)
    export_file_name = muser_data_path(data_file_name)+'.ms'
    # export_blockvisibility_to_ms(export_file_name, vis_list, source_name='SUN')

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

    main(parser.parse_args())