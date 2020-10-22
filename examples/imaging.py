import os
import sys

sys.path.append(os.path.join('..', '..'))

results_dir = '/tmp/'

from matplotlib import pylab

pylab.rcParams['figure.figsize'] = (8.0, 8.0)
pylab.rcParams['image.cmap'] = 'rainbow'

import numpy

from astropy.coordinates import SkyCoord
from astropy import units as u

from matplotlib import pyplot as plt


from rascil.data_models import PolarisationFrame
from astropy.coordinates import EarthLocation, SkyCoord, ITRS
from rascil.processing_components import create_configuration_from_file
from rascil.processing_components.visibility.coalesce import convert_visibility_to_blockvisibility, \
    convert_blockvisibility_to_visibility
from rascil.processing_components.visibility.base import create_blockvisibility_from_ms, create_visibility_from_ms
from rascil.processing_components import create_visibility, show_image, export_image_to_fits, \
    deconvolve_cube, restore_cube, create_named_configuration, create_test_image, \
    create_image_from_visibility, advise_wide_field, invert_2d, predict_2d
from fund import create_configuration
from muser.data_models.parameters import muser_path, muser_data_path, muser_output_path
from rascil.processing_components.visibility.base import export_blockvisibility_to_ms

import logging
from muser.data_models.parameters import muser_path
import os

log = logging.getLogger()
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

pylab.rcParams['figure.figsize'] = (12.0, 12.0)
pylab.rcParams['image.cmap'] = 'rainbow'

location = EarthLocation(lon=115.2505 * u.deg, lat=42.211833333 * u.deg, height=1365.0 * u.m)

# conf_dir = muser_path('configurations')

# antfile = os.path.join(conf_dir, 'muser-1.csv')
# lowcore = create_configuration_from_file(antfile=antfile,
#              mount='altaz', names='MUSER_%d',
#              diameter=2.0, name='MUSER', location=location)

msfile = muser_output_path("ms_flag_test.ms")
# msfile = muser_output_path("CSRH_20151122-125000_100457483.ms")
# msfile = muser_output_path("selfmodel.ms")

ch= numpy.arange(1)
vis = create_blockvisibility_from_ms(msfile) #,  start_chan=52,end_chan=52)

bvis = vis[0]
from rascil.processing_components.flagging.operations import flagging_blockvisibility,flagging_blockvisibility_with_bl
# bvis = flagging_blockvisibility(vis[0], antenna=[8,9,10,11,27])
# baseline = [ [4,0],[4,1],[5,4],[21,4],[24,4],[25,4],[26,4],[27,4],[28,4],[29,4],[30,4],[31,4],[32,4],[36,4],[38,4],[39,4]]
# baseline.append([[17,4],[17,13],[19,17],[26,17],[27,17],[28,17],[29,17],[30,17],[31,17],[39,17]])
# bvis = flagging_blockvisibility_with_bl(vis[0],baseline)
#
# vis_list = []
# vis_list.append(bvis)
# # Output results
# export_file_name = muser_output_path('ms_flag_test') + '.ms'
#
# export_blockvisibility_to_ms(export_file_name, vis_list, source_name='SUN')


# Flagging
vt = convert_blockvisibility_to_visibility(bvis)


advice = advise_wide_field(vt, guard_band_image=3.0, delA=0.1,
                           oversampling_synthesised_beam=4.0)
cellsize = advice['cellsize']

plt.clf()
plt.plot(vt.data['uvw'][:,0], vt.data['uvw'][:,1], '.', color='b')
plt.plot(-vt.data['uvw'][:,0], -vt.data['uvw'][:,1], '.', color='b')
# plt.xlim([-400.0, 400.0])
# plt.ylim([-400.0, 400.0])
plt.show()

# To check that we got the prediction right, plot the amplitude of the visibility.
uvdist=numpy.sqrt(vt.data['uvw'][:,0]**2+vt.data['uvw'][:,1]**2)
plt.clf()
plt.plot(uvdist, numpy.abs(vt.data['vis']), '.')
plt.xlabel('uvdist')
plt.ylabel('Amp Visibility')
plt.show()

model = create_image_from_visibility(vt, cellsize=cellsize, npixel=1024, polarisation_frame=PolarisationFrame('stokesIV'))
dirty, sumwt = invert_2d(vt, model, context='2d')
psf, sumwt = invert_2d(vt, model, context='2d', dopsf=True)

dirty = show_image(dirty)
dirty.show()
show_image(psf)
# print("Max, min in dirty image = %.6f, %.6f, sumwt = %f" % (dirty.data.max(), dirty.data.min(), sumwt))

# print("Max, min in PSF         = %.6f, %.6f, sumwt = %f" % (psf.data.max(), psf.data.min(), sumwt))

# export_image_to_fits(dirty, '%s/imaging_dirty.fits'%(results_dir))
# export_image_to_fits(psf, '%s/imaging_psf.fits'%(results_dir))

comp, residual = deconvolve_cube(dirty, psf, niter=100, threshold=0.001, fractional_threshold=0.001,
                                 gain=0.2,algorithm='msclean',scales=[0, 3, 10 ,30, 50,100])

restored = restore_cube(comp, psf, residual)

# Show the results

fig=show_image(comp)
plt.title('Solution')
fig=show_image(residual)
plt.title('Residual')
fig=show_image(restored)
plt.title('Restored')


