"""Count MUSER Frame Number

"""

from muser.data_models.parameters import muser_path, muser_data_path
from muser.components.utils.installation_checks import check_data_directory
from muser.components.io.muser_data_read import MuserDataReader
from muser.data_models.muser_data import MuserData
from astropy.time import Time
import logging
import os

log = logging.getLogger('logger')


def init_logging():
    logging.basicConfig(filename='check_muser_data.log',
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)


def count_frame_number(time_start, time_end):
    muser = MuserData(sub_array=1)
    if not muser.init_data_environment():
        print("No data environment prepared, exit.")
        return -1
    if not muser.search_first_file(frame_time=time_start):
        print("Cannot find observational data or not a MUSER file.")
        return -1
    print("Checking MUSER File Information V20200801")
    print("First Observational Time (utc): {}".format(muser.current_frame_time.isot))
    print("Filename {} is a valid MUSER Data File.".format(muser.current_file_name))
    # Check data
    if not muser.search_frame(time_start):
        print("Cannot locate the specified frame")
        return -1
    count = 0
    while True:
        if muser.read_full_frame():
            count = count + 1
            if muser.current_frame_time > Time(time_end,format='isot'):
                break
            if not muser.read_one_frame():
                break
        # print("Current Observational Time {}".format(muser.current_frame_time.isot))
        # print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser.is_loop_mode else "Non Loop",
        #                                                      muser.frequency))
        # print("Sub Band: {} - Sub Channel {}".format(muser.sub_band, muser.sub_channels))
        else:
            break
    return count


if __name__ == '__main__':
    # init_logging()
    print("Total Full Frames: %d" % (count_frame_number('2015-11-22T12:45:30', '2015-11-22T12:46:10')))
