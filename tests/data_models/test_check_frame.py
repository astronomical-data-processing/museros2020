"""Check MUSER data frame header

"""

"""Check MUSER Raw data Integrity

"""

import unittest

from muser.data_models.parameters import muser_path, muser_data_path
from muser.components.utils.installation_checks import check_data_directory
from muser.components.io.muser_data_read import MuserDataReader
import logging

log = logging.getLogger('logger')

def init_logging():
    logging.basicConfig(filename='check_muser_data.log',
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)

class TestMuserData(unittest.TestCase):

    def get_muser_data_path(self):
        self.data_path = muser_data_path()
        check_data_directory()

    def test_read_data(self):
        data_file_name = 'CSRH_20151122-093500_89058131'
        mdr = MuserDataReader(1, data_file_name)
        if mdr.search_first_file():
            log.info('Cannot find the data file.')
            exit(1)
        mdr.get_file_info()


if __name__ == '__main__':
    init_logging()
    unittest.main()