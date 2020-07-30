""" Data path testing


"""

import unittest

from muser.data_models.parameters import muser_path, muser_data_path
from muser.components.utils.installation_checks import check_data_directory


class TestDataPath(unittest.TestCase):
    def test_muser_data_path(self):
        data_path = muser_data_path('configurations')
        check_data_directory()


if __name__ == '__main__':
    unittest.main()