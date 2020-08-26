"""Check MUSER data frame header

"""

"""Check MUSER Raw data Integrity

"""

import unittest

from muser.data_models.parameters import muser_path, muser_data_path
from muser.components.utils.installation_checks import check_data_directory
from muser.components.io.muser_data_read import MuserDataReader
import logging
import numpy
from muser.data_models.muser_phase import MuserPhase

log = logging.getLogger('logger')


def init_logging():
    logging.basicConfig(filename='test_strip_stop.log',
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)


class TestReadData(unittest.TestCase):

    def get_muser_data_path(self):
        self.data_path = muser_data_path()
        check_data_directory()

    def test_read_data(self):
        data_file_name = muser_data_path('CSRH_20151122-125000_100457483')
        muser = MuserDataReader(1, file_name=data_file_name)
        if not muser.search_first_file():
            log.info('Cannot find the data file.')
            exit(1)
        muser.skip_frames(11)
        if not muser.read_one_frame():
            log.error("Cannot read proper frame")
            exit(1)
        print("Loading Phase Calibration File")
        phase_cal = MuserPhase(muser.sub_array, muser.is_loop_mode, muser.current_frame_time)
        if not phase_cal.load_calibration_data():
            print("Cannot find phase calibration file. ")
            exit(1)
        print("File shape", phase_cal.phase_data.shape)
        phase_cal_file = numpy.arctan2(phase_cal.phase_data.imag,phase_cal.phase_data.real)
        print(muser.current_frame_time)
        if muser.read_data():
            phase = numpy.arctan2(muser.block_data.imag/muser.block_data.real)
            assert(phase[1,0,0] == -phase[0,1,0])
        


if __name__ == '__main__':
    init_logging()
    unittest.main()
