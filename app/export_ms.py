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
    location = EarthLocation(lon=115.2505, lat=42.211833333, height=1365.0)
    if name == 'MUSER-2':
        lowcore = create_configuration_from_file(antfile="muser-2.csv",
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=2.0, name='MUSER', location=location, **kwargs)
    else:
        lowcore = create_configuration_from_file(antfile="muser-1.csv",
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=4.0, name='MUSER', location=location, **kwargs)
    return lowcore

# def read_muser_data(file_name: str=''):

if __name__ == '__main__':

    if len(sys.argv) == 1:
        muser='MUSER-1'
    else:
        muser='MUSER-2'

    data_file_name = 'CSRH_20151122-093500_89058131'
    file_name = muser_data_path(data_file_name)
    muser = MuserData(sub_array = 1, file_name = file_name)
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


    # MUSER-1 400Mhz - 1975 Mhz
    # MUSER-2 2000 Mhz - 15000 Mhz

    # Create configuration
    muser_core = create_configuration(muser)

    frequency = numpy.array([freq * 1e6])
    # Directory of storing result

    channel_bandwidth = numpy.array([25e6])
    reffrequency = numpy.max(frequency)

    # Create Phase Centre
    phasecentre = SkyCoord(ra=+80 * u.deg, dec=41 * u.deg, frame='icrs', equinox='J2000')

    vt = create_visibility(muser_core, times, frequency, channel_bandwidth=channel_bandwidth,
                           weight=1.0, phasecentre=phasecentre,
                           polarisation_frame=PolarisationFrame('stokesI'))

    lowr3 = create_named_configuration('LOWBD2', rmax=750.0)

    times = numpy.zeros([1])
    frequency = numpy.array([1e8])
    channelbandwidth = numpy.array([1e6])
    phasecentre = SkyCoord(ra=+15.0 * u.deg, dec=-45.0 * u.deg, frame='icrs', equinox='J2000')

    bvis = create_blockvisibility(lowr3, times, frequency, phasecentre=phasecentre,
                                  weight=1.0, polarisation_frame=PolarisationFrame('stokesI'),
                                  channel_bandwidth=channelbandwidth)

    vt = convert_blockvisibility_to_visibility(bvis)

    advice = advise_wide_field(vt, guard_band_image=3.0, delA=0.1, facets=1, wprojection_planes=1,
                               oversampling_synthesised_beam=4.0)
    cellsize = advice['cellsize']

    m31image = create_test_image(frequency=frequency, cellsize=cellsize)
    nchan, npol, ny, nx = m31image.data.shape
    m31image.wcs.wcs.crval[0] = vt.phasecentre.ra.deg
    m31image.wcs.wcs.crval[1] = vt.phasecentre.dec.deg
    m31image.wcs.wcs.crpix[0] = float(nx // 2)
    m31image.wcs.wcs.crpix[1] = float(ny // 2)
    vt = predict_list_serial_workflow([vt], [m31image], context='2d')[0]
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

    v = convert_visibility_to_blockvisibility(vt)
    vis_list = []
    vis_list.append(v)
    export_blockvisibility_to_ms(msoutfile, vis_list, source_name='M31')