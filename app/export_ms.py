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
from rascil.processing_components import weight_visibility
from rascil.processing_components import create_awterm_convolutionfunction
from rascil.processing_components import apply_bounding_box_convolutionfunction

# Use workflows for imaging
from rascil.workflows.rsexecute.execution_support.rsexecute import rsexecute

from rascil.workflows.shared.imaging.imaging_shared import imaging_contexts
from rascil.workflows import predict_list_rsexecute_workflow, \
    invert_list_rsexecute_workflow, deconvolve_list_rsexecute_workflow, \
    residual_list_rsexecute_workflow, restore_list_rsexecute_workflow

import logging
from rascil.data_models.parameters import rascil_path



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
