"""Check MUSER Raw data Integrity

"""


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


def test_read_data():
    data_file_name = 'CSRH_20151122-093500_89058131'
    mdr = MuserDataReader(1, muser_data_path(data_file_name))
    if mdr.search_first_file():
        log.info('Cannot find the data file.')
        exit(1)
    mdr.get_file_info('2015-11-01T11:34:00','2015-11-01T11:34:00')


if __name__ == '__main__':

    test_read_data()