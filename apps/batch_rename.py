# %matplotlib inline

import os
import sys
import numpy
import copy
import matplotlib
import argparse

# from matplotlib import plt.savefig
from astropy.coordinates import EarthLocation, SkyCoord, ITRS, AltAz

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# sys.path.append(os.path.join('..','..'))

from muser.data_models.muser_data import MuserData
from muser.data_models.muser_phase import MuserPhase
from muser.data_models.parameters import muser_path, muser_data_path, muser_output_path
from muser.data_models.parameters import muser_data_list

import logging

def main(args):
    if args.muser == '1':
        muser_array = 'MUSER1'
    else:
        muser_array = 'MUSER2'

    data_path = muser_data_path()
    file_list = muser_data_list(data_path)
    if len(file_list) == 0:
        print("Cannot find any file.")
        return False

    for file_name in file_list:
        full_file_name = muser_data_path(file_name)
        muser = MuserData(sub_array=args.muser, file_name=full_file_name)
        if not muser.open_data_file():
            print("Cannot find observational data or not a MUSER file.")
            continue
        if not muser.read_one_frame():
            print("Cannot read frame information.")
            continue
        current_frame_time = muser.current_frame_time
        muser.close_file()

        new_file_name = ('CSRH_%04d%02d%02d-%02d%02d') % (current_frame_time.datetime.year,
                                                   current_frame_time.datetime.month,
                                                   current_frame_time.datetime.day,
                                                   current_frame_time.datetime.hour,
                                                   current_frame_time.datetime.minute)
        # File name has been modified.
        if new_file_name in file_name:
            continue

        # Rename file name
        new_file_name = muser_data_path(new_file_name)
        print("Rename file %s to %s" % (full_file_name,new_file_name))
        os.rename(full_file_name,new_file_name)


    print("Done. ")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Muser Data Batch Rename')
    parser.add_argument('-m', "--muser", type=int, default=1, help='The MUSER array')
    main(parser.parse_args())
