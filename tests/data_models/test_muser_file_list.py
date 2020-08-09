""" Data file list testing


"""

import unittest

from muser.data_models.parameters import muser_path, muser_data_path, muser_data_list
from muser.components.utils.installation_checks import check_data_directory


class TestFileList(unittest.TestCase):
    def test_muser_data_list(self):
        data_path = muser_data_path()
        file_list = muser_data_list(data_path)
        print(file_list)


if __name__ == '__main__':
    unittest.main()