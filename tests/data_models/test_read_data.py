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
        data_file_name = muser_data_path('CSRH_20151122-094000_89058131') #CSRH_20151122-125000_100457483')
        mdr = MuserDataReader(1, file_name=data_file_name)
        if not mdr.search_first_file():
            log.info('Cannot find the data file.')
            exit(1)
        mdr.skip_frames(2)
        if not mdr.read_one_frame():
            log.error("Cannot read proper frame")
            exit(1)
        print(mdr.current_frame_time)
        if mdr.read_data():
            # mdr.delay_process('satellite')
            phase = numpy.arctan(mdr.block_data.imag/mdr.block_data.real)
            assert(phase[1,0,0] == -2.51372)
            


if __name__ == '__main__':
    init_logging()
    unittest.main()
