''' Phase Calibration

'''

__all__ = ["MuserPhase"]

import logging
from muser.data_models.muser_frame_models import MuserBase, MuserFrame
from muser.data_models.parameters import muser_path, muser_calibration_path, muser_output_path, muser_data_path
from muser.data_models.muser_data import MuserData
from astropy.time import Time
import numpy

log = logging.getLogger('muser')


class MuserPhase(object):
    def __init__(self, sub_array=1, loop_mode=True, obs_date: Time = Time('2015-01-01T00:00:00', format='isot')):
        self.sub_array = sub_array
        self.obs_date = obs_date
        self.is_loop_mode = loop_mode
        if self.sub_array == 1:
            self.antennas = 40
        else:
            self.antennas = 60
        muser_para = {0: (16, 40, 44, 4, 2), 1: (16, 60, 64, 33, 2)}
        self.sub_channels, self.antennas, self.dr_output_antennas, self.frame_number, self.polarization_number = \
            muser_para[self.sub_array - 1]
        self.sub_channels = 16
        if loop_mode:
            self.phase_data = numpy.zeros((1, self.antennas, self.antennas, self.sub_channels * self.frame_number, 2),
                                          dtype='complex')
        else:
            self.phase_data = numpy.zeros((1, self.antennas, self.antennas, self.sub_channels, 1), dtype='complex')

        # self.phase_data = numpy.zeros((4, 2, self.antennas * (self.antennas - 1) // 2, 16), dtype='complex')

    def load_calibration_data(self, file_name=None):

        # if self.last_calibration_priority == self.calibration_priority and self.last_date_time.day == self.start_date_time.day and self.last_date_time.month == self.start_date_time.month and self.last_date_time.year == self.start_date_time.year and self.last_sub_array == self.sub_array:
        #     return
        if file_name is None:
            year = self.obs_date.datetime.year
            month = self.obs_date.datetime.month
            day = self.obs_date.datetime.day
            file_name = muser_calibration_path("MUSER%1d-%04d%02d%02d.CAL" % (self.sub_array, year, month, day))

        import os
        if os.path.isfile(file_name):
            caldata = numpy.fromfile(file_name, dtype=complex)
            if self.is_loop_mode:
                self.phase_data = caldata.reshape(1, self.antennas, self.antennas, self.sub_channels * self.frame_number, 2)
            else:
                self.phase_data = caldata.reshape(1, self.antennas, self.antennas, self.sub_channels, 1)
            log.debug("Load Calibrated data.")
            return True
        else:
            log.error("Cannot find calibrated data.")
            return False
