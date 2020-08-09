"""Export 

"""

import os, sys

path = os.path.abspath(os.path.dirname(__file__))

from argparse import *
import datetime
import numpy
import traceback
import logging
from muser.data_models.muser_data import MuserData
from muser.data_models.parameters import muser_path, muser_data_path
from rascil.processing_components.visibility.coalesce import convert_visibility_to_blockvisibility, \
    convert_blockvisibility_to_visibility
from astropy.time import Time

log = logging.getLogger('muser')

class Phase:
    def __init__(self, sub_array=None, is_loop_mode=None, obs_date_time=None, frame_number=None, file_name=None):
        self.sub_array = sub_array
        self.frame_number = frame_number
        self.data_source = 0
        if obs_date_time is not None:
            self.obs_date_time = Time(obs_date_time, format='isot')
        else:
            self.obs_date_time = None
        self.is_loop_mode = is_loop_mode

        if obs_date_time is not None and len(obs_date_time)>0:
            log.info("Searched phase data date and time: %s" % obs_date_time)

        self.file_name = file_name.strip()

    def calibration(self):
        if self.frame_number < 0:
            log.error("You should input a positive number!")
            return False
        file_name = muser_data_path(self.file_name)
        muser_calibration = MuserData(self.sub_array, file_name)
        if self.obs_date_time is not None:
            muser_calibration.set_data_date_time(self.obs_date_time)
        log.info('Reading Visibility Data of calibration......')
        if not muser_calibration.open_data_file():
            print("Cannot find observational data or not a MUSER file.")
            exit(1)

        muser_calibration.skip_frames(self.frame_number)
        self.last_sub_band = -1
        self.last_polarization = -1

        if self.is_loop_mode == True:
            frame_NUM = muser_calibration.frame_number * 2

            calibration_Data = numpy.ndarray(
                shape=(muser_calibration.frame_number, muser_calibration.polarization_number,
                       muser_calibration.antennas * (muser_calibration.antennas - 1) // 2, 16),
                dtype=complex)
        else:
            frame_NUM = 1
            calibration_Data = numpy.ndarray(
                shape=(muser_calibration.antennas * (muser_calibration.antennas - 1) // 2, 16),
                dtype=complex)

        for i in range(frame_NUM):
            if (muser_calibration.read_one_frame() == False):  # 32*8bits
                log.error('Cannot read a frame.')
                return False
            muser_calibration.read_data()

            self.year = muser_calibration.current_frame_time.datetime.year
            self.month = muser_calibration.current_frame_time.datetime.month
            self.day = muser_calibration.current_frame_time.datetime.day

            log.info("Reading No. %d %s %d %d" % (i, muser_calibration.current_frame_time.isot, muser_calibration.sub_band, muser_calibration.polarization))
            #According to the DISCUSSION with LJ.Chen
            #muser_calibration.delay_process('satellite')
            # Delay processing for satellite
            if self.sub_array == 2:
                if muser_calibration.current_frame_header.strip_switch == 0xCCCCCCCC:
                    muser_calibration.delay_process("satellite")

            self.last_sub_band =  muser_calibration.sub_band
            self.last_polarization = muser_calibration.polarization

            bl = 0
            for antenna1 in range(0, muser_calibration.antennas - 1):
                for antenna2 in range(antenna1 + 1, muser_calibration.antennas):
                    for channel in range(0, muser_calibration.sub_channels):
                        if self.is_loop_mode == True:
                            calibration_Data[muser_calibration.sub_band][muser_calibration.polarization][bl][channel] = muser_calibration.baseline_data[bl][channel]
                        else:
                            calibration_Data[bl][channel] = muser_calibration.baseline_data[bl][channel]
                    bl = bl + 1

        file_name = muser_data_path("MUSER%1d-%04d%02d%02d.CAL"%(self.sub_array,self.year, self.month, self.day))
        log.info("Writing to file: " + os.path.basename(file_name))
        calibration_Data.tofile(file_name)
        log.info("Export phase done.")
        return True

def export_phase(args):
    muser = args.muser
    input_file = args.file
    frame = args.number
    start = args.start
    loop_mode = args.loop

    cal = Phase(muser, loop_mode, start, frame, input_file)
    cal.calibration()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='List Muser Data Information for Each Frame')
    parser.add_argument('-m', "--muser", type=int, default=1, help='The MUSER array')
    parser.add_argument('-f', "--file", type=str, required=True, default='', help='The file name')
    parser.add_argument('-l', "--loop", type=bool, default=True, help='Loop Mode')
    parser.add_argument('-s', "--start", type=str, default=None, help='The beginning time')
    parser.add_argument('-c', "--calib", type=str, default='', help='The beginning time')
    parser.add_argument('-n', "--number", type=int, default=1, help='The number of frames')

    export_phase(parser.parse_args())





