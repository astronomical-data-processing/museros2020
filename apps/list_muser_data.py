#! /usr/bin/env python

"""List Muser Data Frame Information

"""
from muser.data_models.parameters import muser_path, muser_data_path
from muser.components.utils.installation_checks import check_data_directory
from muser.components.io.muser_data_read import MuserDataReader
from muser.data_models.muser_data import MuserData
import logging
import os

log = logging.getLogger('logger')

def list_muser_data(args):
    loop_number = args.line
    data_file_name = args.file
    # data_file_name = 'CSRH_20151122-093500_89058131'
    file_name = muser_data_path(data_file_name)
    muser = MuserData(sub_array = 1, file_name = file_name)
    if not muser.open_data_file():
        print("Cannot find observational data or not a MUSER file.")
        exit(1)
    # Read first Frame
    muser.read_one_frame()
    # Locate a specified frame
    print("Checking MUSER File Information V20200801")
    print("First Observational Time {}".format(muser.current_frame_time.isot))
    print("Filename {} is a valid MUSER Data File.".format(file_name))
    print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser.is_loop_mode else "Non Loop", muser.frequency))
    # Check data
    for i in range(loop_number):
        print("No.{} Time:{} Pol:{} SubBand:{} Freq:{}".format(i,muser.current_frame_time.isot,muser.polarization,muser.sub_band,muser.frequency))
        muser.read_one_frame()
    # print("Sub Band: {} - Sub Channel {}".format(muser.sub_band,muser.sub_channels))
    muser.close_file()

if __name__ == '__main__':
    # Usage: -l 10 -f CSRH_20151122-093500_89058131

    import argparse

    parser = argparse.ArgumentParser(description='List Muser Data Information for Each Frame')
    parser.add_argument('-f', "--file", type=str, default='', help='The file name')
    parser.add_argument('-l', "--line", type=int, default=1, help='The number of frames')

    list_muser_data(parser.parse_args())