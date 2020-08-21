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
from muser.data_models.parameters import muser_path, muser_data_path, muser_calibration_path, muser_output_path
from rascil.processing_components.visibility.coalesce import convert_visibility_to_blockvisibility, \
    convert_blockvisibility_to_visibility
from astropy.time import Time

log = logging.getLogger('muser')


class Phase:
    def __init__(self, sub_array=None, start_time=None, end_time=None, output_file=None, fringe=False):
        self.sub_array = sub_array
        self.data_source = 0
        if start_time is not None:
            self.start_time = start_time
        else:
            self.start_time = None
        if end_time is not None:
            self.end_time = end_time
        else:
            self.end_time = None
        if output_file is not None:
            self.output_file = output_file
        else:
            self.output_file = None
        self.fringe_stop = fringe

    def calibration(self):
        muser_calibration = MuserData(sub_array=self.sub_array)
        log.info('Reading Visibility Data of calibration......')
        if not muser_calibration.init_data_environment():
            print("No data environment prepared, exit.")
            return -1
        if not muser_calibration.search_first_file(frame_time=self.start_time):
            print("Cannot find observational data or not a MUSER file.")
            return -1

        data_file_name = muser_calibration.current_file_name
        print("Checking MUSER File Information V20200801")
        print("First Observational Time {}".format(muser_calibration.current_frame_time.isot))
        # Check data
        print("Filename {} is a valid MUSER Data File.".format(data_file_name))
        print("Current Observational Time UTC:  {}".format(muser_calibration.current_frame_utc_time.isot))
        print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser_calibration.is_loop_mode else "Non Loop",
                                                             muser_calibration.frequency))
        print("First frame Sub Band: {} - Sub Channel {}".format(muser_calibration.sub_band, muser_calibration.sub_channels))

        # count total frames
        muser_calibration.search_frame(search_time=self.start_time)

        self.last_sub_band = -1
        self.last_polarization = -1

        if muser_calibration.is_loop_mode == True:
            if self.end_time is None:
                total_frames = 1
            else:
                total_frames = muser_calibration.count_frame_number(self.start_time, self.end_time)
            self.block_full_data = numpy.zeros(
                [muser_calibration.antennas, muser_calibration.antennas,
                 muser_calibration.sub_channels * muser_calibration.frame_number,
                 2], dtype='complex')
        else:
            if self.end_time is None:
                total_frames = 1
            else:
                total_frames = muser_calibration.count_frame_number(self.start_time, self.end_time)
            self.block_full_data = numpy.zeros(
                [muser_calibration.antennas, muser_calibration.antennas, muser_calibration.channels],
                dtype='complex')

        # Re-Search file
        # if not muser_calibration.search_first_file(frame_time=self.start_time):
        #     print("Cannot find observational data or not a MUSER file.")
        #     return -1
        if not muser_calibration.search_frame(self.start_time):
            print("Cannot locate the specified frame")
            return -1

        self.year = muser_calibration.current_frame_time.datetime.year
        self.month = muser_calibration.current_frame_time.datetime.month
        self.day = muser_calibration.current_frame_time.datetime.day

        count = 0
        while count < total_frames:
            if not muser_calibration.read_full_frame(read_data=True):
                print("File reading error. ")
                exit(1)

            print("Reading No. %d %s" % (
                count, muser_calibration.first_frame_time.isot))
            log.info("Reading No. %d %s" % (
                count, muser_calibration.first_frame_time.isot))

            # Delay processing for satellite
            if self.fringe_stop:
                if self.sub_array == 2:
                    if muser_calibration.current_frame_header.strip_switch == 0xCCCCCCCC:
                        muser_calibration.delay_process("satellite")
                else:
                    muser_calibration.delay_process('satellite')

            self.last_sub_band = muser_calibration.sub_band
            self.last_polarization = muser_calibration.polarization

            # inject data
            if muser_calibration.is_loop_mode:
                self.block_full_data[:, :, :] += muser_calibration.block_full_data[:, :, :]
            else:
                self.block_full_data[:, :, :] += muser_calibration.block_data[:, :, :]
            count = count + 1
            # bl = 0
            # for antenna1 in range(0, muser_calibration.antennas - 1):
            #     for antenna2 in range(antenna1 + 1, muser_calibration.antennas):
            #         for channel in range(0, muser_calibration.sub_channels):
            #             if self.is_loop_mode == True:
            #                 calibration_Data[muser_calibration.sub_band][muser_calibration.polarization][bl][channel] = \
            #                     muser_calibration.baseline_data[bl][channel]
            #             else:
            #                 calibration_Data[bl][channel] = muser_calibration.baseline_data[bl][channel]
            #         bl = bl + 1

        # Mean Value
        self.block_full_data /= count
        if self.output_file is None:
            file_name = muser_calibration_path(
                "MUSER%1d-%04d%02d%02d.CAL" % (self.sub_array, self.year, self.month, self.day))
        else:
            file_name = muser_calibration_path(self.output_file)
        log.info("Writing to file: " + os.path.basename(file_name))
        self.block_full_data.tofile(file_name)
        log.info("Export phase done.")
        return True


def export_phase(args):
    muser = args.muser
    start = args.start
    end_time = args.end
    if len(args.output) != 0:
        output_file = args.output
    else:
        output_file = None
    fringe = args.fringe
    cal = Phase(sub_array=muser, start_time =start, end_time =end_time, output_file =output_file, fringe=fringe)
    cal.calibration()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='List Muser Data Information for Each Frame')
    parser.add_argument('-m', "--muser", type=int, default=1, help='The MUSER array')
    parser.add_argument('-s', "--start", type=str, default=None, help='The beginning time')
    parser.add_argument('-e', "--end", type=str, default=None, help='The end time')
    parser.add_argument('-f', "--fringe", type=bool, default=False, help='Fringe Stop')
    parser.add_argument('-o', "--output", type=str, default='', help='The output filename')

    export_phase(parser.parse_args())
