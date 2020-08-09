"""Check MUSER Raw data Integrity

"""


from muser.data_models.parameters import muser_path, muser_data_path
from muser.components.utils.installation_checks import check_data_directory
from muser.components.io.muser_data_read import MuserDataReader
from muser.data_models.muser_data import MuserData
import logging
import os

log = logging.getLogger('logger')

def init_logging():
    logging.basicConfig(filename='check_muser_data.log',
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)

def check_data_info():
    data_file_name = 'CSRH_20151122-093500_89058131'
    file_name = muser_data_path(data_file_name)
    muser = MuserData(sub_array = 1, file_name = file_name)
    if not muser.search_first_file():
        print("Cannot find observational data or not a MUSER file.")
        exit(1)
    print("Checking MUSER File Information V20200801")
    print("First Observational Time (utc): {}".format(muser.current_frame_time.isot))
    # Check data
    muser.search_frame('2015-11-22T09:41:00')
    print("Filename {} is a valid MUSER Data File.".format(file_name))
    print("Current Observational Time {}".format(muser.current_frame_time.isot))
    print("Observational Mode: {} \nFrequency {}".format("LOOP" if muser.is_loop_mode else "Non Loop", muser.frequency))
    print("Sub Band: {} - Sub Channel {}".format(muser.sub_band,muser.sub_channels))
    muser.read_data()
    print(muser.baseline_data)
    print(muser.block_data)

if __name__ == '__main__':
    # init_logging()
    check_data_info()