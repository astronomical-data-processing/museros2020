""" MUSER Data Reader


"""

__all__ = ['MuserDataReader']

from astropy.coordinates import EarthLocation, SkyCoord

import os

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

from rascil.data_models.parameters import rascil_path

import logging,os
import numpy
from muser.data_models.muser_data import MuserData

from muser.data_models.parameters import muser_path, muser_data_path

from astropy.time import Time

log = logging.getLogger('logger')

class MuserDataReader(MuserData):

    def __init__(self, sub_array=1, file_name = None):
        """
        Main function call. Process raw data: delay process and sum
        """
        super(MuserDataReader, self).__init__(sub_array, file_name=file_name)

        # self.config = os.path.join(os.path.abspath(os.path.dirname(__file__)), "museruvfits.xml")

        # muserRawData Class

    def __del__(self):
        if self.in_file is not None:
            self.close_file()

    def get_data_info(self, start_time, end_time, integral_period, realtime=False):
        s_time = Time(start_time,format='isot', scale='utc')
        self.set_data_date_time(s_time)
        (self.t_len, self.chan_len, self.bl_len, self.pol_len, self.ri_len) = (
            1, self.sub_channels, self.antennas * (self.antennas - 1) / 2, 1, 2)

        # If cannot locate a proper frame, return with False
        if realtime:
            if self.search_frame_realtime() == False:
                log.error("cannot find observational data.")
                return False, None, None, None
        else:
            if self.search_first_frame() == False:
                log.error("cannot find observational data.")
                return False, None, None, None

        frame = [2500, 20625]  # 25,206.25 -> 2500, 20625, so integraltime * 1e8 instead of 1e6

        if self.is_loop_mode == False:
            integral_number = int(integral_period * 1e6 // 3125)  # 3.125 -> 3125
        else:
            if ((integral_period * 1e8) % frame[self.sub_array - 1]) == 0:
                integral_number = integral_period * 1e5 // frame[self.sub_array - 1]
            else:
                integral_number = int(integral_period * 1e5 // frame[self.sub_array - 1]) + 1

        if ((start_time.get_date_time() - end_time.get_date_time()).seconds * 1E5) % (
                    integral_number * (312.5 if self.is_loop_mode == False else frame[self.sub_array - 1])) == 0:
            total_loop_number = int(abs((end_time.get_date_time() - start_time.get_date_time()).seconds * 1E5) / (
                integral_number * (312.5 if self.is_loop_mode == False else frame[self.sub_array - 1])))
        else:
            total_loop_number = int(
                abs((end_time.get_date_time() - start_time.get_date_time()).seconds * 1E5) / (
                    integral_number * (312.5 if self.is_loop_mode == False else frame[self.sub_array - 1]))) + 1

        return self.current_frame_time.get_date_time(), integral_number, total_loop_number, self.is_loop_mode


    def get_frame_info(self, start_time, end_time):
        s_time = Time(start_time, format='isot', scale='utc')
        self.set_data_date_time(start_time)
        (self.t_len, self.chan_len, self.bl_len, self.pol_len, self.ri_len) = (
            1, self.sub_channels, self.antennas * (self.antennas - 1) / 2, 1, 2)

        # If cannot locate a proper frame, return with False
        if self.search_first_frame() == False:
            log.error("cannot find observational data.")
            return False, None, None, None

        frame = [2500, 20625]  # 25,206.25 -> 2500, 20625, so integraltime * 1e8 instead of 1e6
        integral_number = 1

        if ((start_time.get_date_time() - end_time.get_date_time()).seconds * 1E5) % (
                    integral_number * (312.5 if self.is_loop_mode == False else frame[self.sub_array - 1])) == 0:
            total_loop_number = int(abs((end_time.get_date_time() - start_time.get_date_time()).seconds * 1E5) / (
                integral_number * (312.5 if self.is_loop_mode == False else frame[self.sub_array - 1])))
        else:
            total_loop_number = int(
                abs((end_time.get_date_time() - start_time.get_date_time()).seconds * 1E5) / (
                    integral_number * (312.5 if self.is_loop_mode == False else frame[self.sub_array - 1]))) + 1

        return total_loop_number



    def get_file_info(self, start_time, end_time): # 2015-11-01 11:34:00  2015-11-01 11:36:00

        file_info=[]
        s_time = Time(start_time, format='isot', scale='utc')
        self.set_data_date_time(s_time)

        if self.search_first_file() == False:
            log.info("Cannot find observational data.")
            return

        if self.open_raw_file(self.current_file_name) == False:
            log.info("Cannot open observational data.")
            return

        if self.read_one_frame() == False:
            log.info("Error reading frame.")
            return

        file_info.append(self.current_file_name)

        while True:

            if self.open_next_file(1) == False:
                break
            self.read_one_frame()
            if self.current_frame_time.get_time_stamp() > end_time.get_time_stamp():
                break
            file_info.append(self.current_file_name)

        return file_info

    def get_visdata(self, filename, offset, repeat_number):

        self.in_file = open(filename, 'rb')
        self.in_file.seek(offset, 0)

        vis_file= ""
        uvw_file= ""
        date_file = ""
        uvws = []
        visdata = []
        date_FILE = numpy.ndarray(shape = (repeat_number, 7), dtype= float)

        for iLoop in range(0, repeat_number):

            # Read Observational Data
            self.read_one_frame()
            self.read_data()
            if self.debug:
                log.info("Time: %s BAND:%5d POL: %-3s" % (
                self.current_frame_date_time, self.channel_group, "LL" if self.polarization == 1 else "RR"))

            # Skip the first frame which is right polarization in loop mode

            # Delay process and Strip stop
            if self.current_frame_header.strip_switch == 0xCCCCCCCC:
                if self.debug:
                    log.debug("Strip rotation correction")
                self.delay_process('sun')

            if self.no_calibration == 1:
                # Load Calibration file from disk according to the observational date and time
                self.load_calibration_data()
                # Calibration
                self.calibration()

            visdata.append(self.baseline_data)

            if iLoop == 0:
                file_VIS = ('%4d%02d%02d-%02d%02d%02d_%03d%03d%03d.vis') % (   #FILENAME VIS UVW
                    self.current_frame_time.get_detail_time())
                vis_file = self.env.vis_file(self.sub_array, file_VIS)

                file_UVW = ('%4d%02d%02d-%02d%02d%02d_%03d%03d%03d.uvw') % (
                    self.current_frame_time.get_detail_time())
                uvw_file = self.env.uvw_file(self.sub_array, file_UVW)

                file_DATE = ('%4d%02d%02d-%02d%02d%02d_%03d%03d%03d.date') % (
                    self.current_frame_time.get_detail_time())
                date_file = self.env.uvw_file(self.sub_array, file_DATE)

            # Sun = self.pyuvfits.makeSource(name=self.obs_target)
            # self.source = Sun
            #
            # obs = Observatory(lon=self.longitude, lat=self.latitude, altitude=self.altitude)
            #
            # array_geometry = self.pyuvfits.ant_array()
            # antenna_array = Array(lat=self.latitude, long=self.longitude, elev=self.altitude, antennas=array_geometry)
            # self.source.midnightJD, midnightMJD = self.pyuvfits.ephem.convert_date(self.obs_date, '00:00:00')
            # self.source.compute(cobs=obs, cdate=self.obs_date, ctime=self.obs_time)


            self.baseline = []
            bl_len = int(self.antennas * (self.antennas - 1) / 2)
            (bl_order, baselines) = self.pyuvfits.config_baseline_ID(bl_len)

            for baseline in baselines:
                vector = baseline[1]
                self.baseline.append(baseline[0])
                if self.hourangle==999 and self.declination ==999:
                    H, d = (self.source.gast - self.source.appra, self.source.appdec)
                else:
                    H, d = self.hourangle, self.declination

            if self.is_loop_mode == True:
                if iLoop %2 == 0:
                    uvws.append(self.pyuvfits.computeUVW(vector, H * 15., d) / light_speed )  # units: SECOND
                    self.obs_date_sum = self.source.midnightJD
                    self.obs_time_sum = self.source.JD - self.source.midnightJD
                    date_FILE[iLoop][0] = self.obs_date_sum
                    date_FILE[iLoop][1] = self.obs_time_sum
                    date_FILE[iLoop][2] = self.source.appra
                    date_FILE[iLoop][3] = self.source.appdec
                    date_FILE[iLoop][4] = self.source.topora
                    date_FILE[iLoop][5] = self.source.topodec
                    date_FILE[iLoop][6] = self.frequency

            else:
                light_speed = numpy.c
                uvws.append(self.pyuvfits.computeUVW(vector, H * 15., d) / light_speed )  # units: SECOND
                # self.obs_date_sum = self.source.midnightJD
                # self.obs_time_sum = self.source.JD - self.source.midnightJD
                # date_FILE[iLoop][0] = self.obs_date_sum
                # date_FILE[iLoop][1] = self.obs_time_sum

        UVW= numpy.array(uvws)
        VIS_DATA = numpy.array(visdata)
        UVW.tofile(vis_file)
        VIS_DATA.tofile(uvw_file)
        date_FILE.tofile(date_file)

        log.info("Writing to file: ", os.path.basename(vis_file),os.path.basename(uvw_file),os.path.basename(date_file))
        return [vis_file, uvw_file, date_file]







