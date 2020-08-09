"""

"""
__all__ = ['MuserData']

import sys, os, datetime, time, math
from muser.data_models.muser_frame_models import MuserBase, MuserFrame
from muser.data_models.parameters import muser_data_list, IF_BANDWIDTH, LOOP_MODE_LOW, NON_LOOP_MODE_LOW, \
    LOOP_MODE_HIGH
from astropy.time import Time
from astropy import units as u
import numpy
import logging
from muser.data_models.parameters import muser_path, muser_data_path

log = logging.getLogger('muser')


class MuserData(MuserFrame):
    def __init__(self, sub_array=1, file_name = None):

        super(MuserData, self).__init__(sub_array=1)
        if file_name is not None:
            self.input_file_name = file_name
        self.if_time_setup = False
        self.if_file_opened = False
        self.in_file = None
        self.integral_time = 0

        self.last_calibration_priority = -1
        self.last_date_time = 0
        self.last_sub_array = -1

        self.current_file_name = ''
        self.file_list = []
        # self.start_date_time = MuserTime()
        # self.initial_frame_time = MuserTime()

    def __del__(self):
        self.close_raw_file()

    # @property
    # def filename(self, name):
    #     self.filename = name

    def set_data_date_time(self, time):
        self.start_date_time = time
        self.initial_frame_time = time

        # set True to if_time_setup
        self.if_time_setup = True

    def set_priority(self, priority):
        self.calibration_priority = priority

    def init_data_environment(self):
        data_path = muser_data_path()
        self.file_list = muser_data_list(data_path)
        if len(self.file_list) == 0:
            return False
        return True

    def get_info(self):
        """Get Muser Data File Information

        """
        if self.filename is None or len(self.filename.strip()) == 0:
            return False
        if self.in_file.closed:
            if self.open_data_file() == False:
                return False

    def search_first_file(self, frame_time = ''):
        '''
        search the first file of processing
        '''
        if len(frame_time) == 0:
            return False
        self.initial_frame_time = Time(frame_time,format='isot')

        if self.open_data_file() == False:
            log.error("Cannot open observational data.")
            return False

        # Read first frame and check the observational date and tim  e
        if self.read_one_frame() == False:
            log.error("Cannot read a frame.")
            return False

        if self.current_file_name != '':  # Start!!!
            self.start_date_time = self.current_frame_time

        # The date and time in the first frame should be
        if self.current_frame_time.datetime.second == 0 and (self.current_frame_time.datetime.microsecond < 312500):
            return True

        if self.start_date_time < self.current_frame_time:
            # previous 1 minute
            log.debug("Change observational file.")
            if self.open_next_file(-1) == False:
                return False
            # Read first frame and check the observational date and time
            if self.read_one_frame() == False:
                self.in_file.close()
                return False
            log.debug("Find data file: %s" % (os.path.basename(self.current_file_name)))
            return False
        else:
            log.debug("Found observational data file: %s" % (os.path.basename(self.current_file_name)))
            return True

    def open_next_file(self, time_minute=1):
        search_date_time = self.initial_frame_time + time_minute * u.minute

        full_file_name = self.muser_data_file_name(self.search_date_time.datetime.year,
                                                   self.search_date_time.datetime.month,
                                                   self.search_date_time.datetime.day,
                                                   self.search_date_time.datetime.hour,
                                                   self.search_date_time.datetime.minute)
        for file in self.file_list:
            if full_file_name in file:
                return self.open_raw_file(file)
        return False

    def muser_data_file_name(self, year, month, day, hour, minute):
        file_name = ('CSRH_%04d%02d%02d-%02d%02d') % (year, month, day, hour, minute)
        for file in self.file_list:
            if file_name in file:
                return muser_data_path(file)
        return ''

    def check_next_file(self):
        '''
        Check file to determine whether the file should be changed
        '''
        offset = [100000, 204800]
        pos = self.in_file.tell()
        if ((pos // offset[self.sub_array - 1]) + 1) >= 19200:
            return True
        else:
            return False

    def open_raw_file(self, file_name):
        if not os.path.exists(file_name):
            log.error("Cannot find file: %s " % (file_name))
            return False
        try:
            self.in_file = open(file_name, 'rb')
            self.in_file.seek(0, 0)
            self.if_read_first_frame_time = False
            self.current_file_name = file_name
            log.debug("File opened: %s" % (os.path.basename(file_name)))
        except:
            # self.in_file.close()
            return False
        return True

    def close_file(self):
        self.close_raw_file()

    def close_raw_file(self):
        try:
            self.in_file.close()
        finally:
            pass

    def open_data_file(self, file_name=None):
        if self.input_file_name == '':
            full_file_name = self.muser_data_file_name(self.initial_frame_time.datetime.year, self.initial_frame_time.datetime.month,
                                                  self.initial_frame_time.datetime.day, self.initial_frame_time.datetime.hour,
                                                  self.initial_frame_time.datetime.minute)
        else:
            full_file_name = self.input_file_name

        return self.open_raw_file(full_file_name)

    def search_first_frame(self):
        '''
        search first frame with proper date and time
        '''

        if self.search_first_file() == False:
            log.error("Cannot find a proper frame from the observational data.")
            return False

        frame_date_time = self.current_frame_time
        if self.start_date_time < self.current_frame_time:
            t_offset = self.current_frame_time - self.start_date_time
        else:
            t_offset = self.start_date_time - self.current_frame_time

        log.debug('Current frame time: %s' % (self.current_frame_time.isot))
        time_offset = t_offset.to_value('s') * 1e6

        skip_frame_number = int(time_offset / 3125) - 1

        log.debug('Time interval %d, skip frames: %d' % (time_offset, skip_frame_number))

        if (skip_frame_number >= 2):
            self.skip_frames(skip_frame_number)

        # If can find, search first frame
        while True:
            if self.is_loop_mode == False:
                if self.start_date_time <= self.current_frame_time:
                    break
            else:
                if self.is_loop_mode == True and self.start_date_time <= self.current_frame_time and self.sub_band == 0 and self.polarization == 0:  # Find file in previous 1 minute
                    break
                if self.is_loop_mode == False and self.start_date_time <= self.current_frame_time:
                    break

            if self.read_one_frame() == False:
                self.in_file.close()
                return False
        log.debug('Frame located.')
        return True

    def search_frame_realtime(self, specified_file=False):
        '''
        search first frame with proper date and time
        '''
        # Search first possible file

        if self.search_first_file() == False:
            log.error("Cannot find a proper frame from the observational data.")
            return False

        frame_date_time = self.current_frame_time
        if self.start_date_time < self.current_frame_time:
            t_offset = self.current_frame_time - self.start_date_time
        else:
            t_offset = self.start_date_time - self.current_frame_time

        # Estimate the number of skip
        # print('current frame time: %04d-%02d-%02d %02d:%02d:%02d %03d%03d%03d' % (
        #     self.current_frame_time.get_detail_time()))
        log.debug('Current frame time: %s' % (self.current_frame_time.isot))
        time_offset = t_offset.datetime.second * 1e6 + t_offset.datetime.microsecond

        skip_frame_number = int(time_offset / 3125) - 1
        log.debug('Time interval %d, skip frames: %d' % (time_offset, skip_frame_number))

        if (skip_frame_number >= 2):
            self.skip_frames(skip_frame_number)

        # If can find, search first frame
        while True:
            if self.is_loop_mode == False:
                if self.start_date_time <= self.current_frame_time:
                    break
            else:

                if specified_file == False:
                    if self.is_loop_mode == False and self.start_date_time <= self.current_frame_time:
                        break
                else:
                    if self.is_loop_mode == True and self.start_date_time <= self.current_frame_time:  # Find file in previous 1 minute
                        break

            if self.read_one_frame() == False:
                self.in_file.close()
                return False
        log.debug('Frame located.')
        return True

    def delay_process(self, planet):
        # print "delay processing..."
        parameter = 0.
        delay = numpy.ndarray(shape=(self.dr_output_antennas), dtype=float)

        if self.sub_array == 1:  # muser-1
            if planet == 'sun':
                parameter = 12.5
            elif planet == 'satellite':
                parameter = 2.5
            delayns = self.delay_compensation.get_delay_value(self.sub_array,
                                                              self.current_frame_time)
            delay = self.par_delay * (10 ** 9) - delayns
        else:  # muser-2
            parameter = 12.5
            delay = self.par_delay

        for channel in range(0, self.sub_channels):
            bl = 0
            for antenna1 in range(0, self.antennas - 1):  # SubChannelsLow = 16
                for antenna2 in range(antenna1 + 1, self.antennas):
                    tg = delay[antenna2] - delay[antenna1]
                    tg0 = int(delay[antenna2]) - int(delay[antenna1])
                    if self.sub_array == 1:
                        Frf = (self.frequency * 1e-6 + channel * 25 + parameter) / 1000.0
                        Fif = (channel * 25 + parameter + 50.0) / 1000.0
                        phai = 2 * numpy.pi * (Frf * tg - Fif * tg0)
                        self.block_data[0, antenna1, antenna2, channel, 0] = complex(
                            self.block_data[0, antenna1, antenna2, channel, 0].real * numpy.cos(phai) +
                            self.block_data[0, antenna1, antenna2, channel, 0].imag * numpy.sin(phai),
                            self.block_data[0, antenna1, antenna2, channel, 0].imag * numpy.cos(phai) -
                            self.block_data[0, antenna1, antenna2, channel, 0].real * numpy.sin(phai))
                        self.baseline_data[bl][channel] = complex(
                            self.baseline_data[bl][channel].real * numpy.cos(phai) +
                            self.baseline_data[bl][channel].imag * numpy.sin(phai),
                            self.baseline_data[bl][channel].imag * numpy.cos(phai) -
                            self.baseline_data[bl][channel].real * numpy.sin(phai))
                    else:
                        Frf = (self.frequency * 1e-6 + (15 - channel) * 25 + parameter) / 1000.0
                        Fif = (channel * 25 + parameter + 50.0) / 1000.0  # local frequency(GHz)
                        phai = 2 * numpy.pi * (-Frf * tg - Fif * tg0)
                        # phai = 2 * pi * Fif * tg0 + 2 * pi * Frf * (tg - tg0)
                        self.block_data[0, antenna1, antenna2, channel, 0] = complex(
                            self.block_data[0, antenna1, antenna2, channel, 0].real * numpy.cos(phai) +
                            self.block_data[0, antenna1, antenna2, channel, 0].imag * numpy.sin(phai),
                            self.block_data[0, antenna1, antenna2, channel, 0].imag * (-1) * numpy.cos(phai) +
                            self.block_data[0, antenna1, antenna2, channel, 0].real * numpy.sin(phai))

                        self.baseline_data[bl][channel] = complex(
                            self.baseline_data[bl][channel].real * numpy.cos(phai) +
                            self.baseline_data[bl][channel].imag * numpy.sin(phai),
                            self.baseline_data[bl][channel].imag * (-1) * numpy.cos(phai) +
                            self.baseline_data[bl][channel].real * numpy.sin(phai))
                    bl = bl + 1

        log.debug("Delay Process and fringe stopping... Done.")

    def phase_calibration(self, cal):
        log.debug("Satellite phase correction")

        if self.sub_array == 1:
            for chan in range(0, self.sub_channels):
                bl = 0
                for antenna1 in range(0, self.antennas - 1):
                    for antenna2 in range(antenna1 + 1, self.antennas):
                        A = numpy.sqrt(
                            self.baseline_data[bl][chan].imag * self.baseline_data[bl][chan].imag +
                            self.baseline_data[bl][chan].real * self.baseline_data[bl][chan].real)

                        phai_sun = numpy.arctan2(self.baseline_data[bl][chan].imag,
                                                 self.baseline_data[bl][chan].real)
                        if self.is_loop_mode == True:
                            phai = phai_sun - numpy.arctan2(
                                cal[self.sub_band][self.polarization][bl][chan].imag,
                                cal[self.sub_band][self.polarization][bl][chan].real)
                        else:
                            phai = phai_sun - numpy.arctan2(cal[bl][chan].imag, cal[bl][chan].real)
                        self.baseline_data[bl][chan] = complex(A * numpy.cos(phai), A * numpy.sin(phai))
                        self.block_data[0, antenna1, antenna2, chan, 0] = self.baseline_data[bl][chan]
                        bl = bl + 1
        else:

            for chan in range(0, self.sub_channels):
                bl = 0
                for antenna1 in range(0, self.antennas - 1):
                    for antenna2 in range(antenna1 + 1, self.antennas):
                        A = numpy.sqrt(
                            self.baseline_data[bl][chan].imag * self.baseline_data[bl][chan].imag +
                            self.baseline_data[bl][chan].real * self.baseline_data[bl][chan].real)
                        phai_sun = numpy.arctan2(self.baseline_data[bl][chan].imag, self.baseline_data[bl][chan].real)
                        if self.is_loop_mode:
                            phai = phai_sun - numpy.arctan2(
                                self.cal[self.sub_band][self.polarization][bl][chan].imag,
                                self.cal[self.sub_band][self.polarization][bl][chan].real)
                        else:
                            phai = phai_sun - numpy.arctan2(self.cal[bl][chan].imag, self.cal[bl][chan].real)
                        self.baseline_data[bl][chan] = complex(A * numpy.cos(phai), A * numpy.sin(phai))
                        self.block_data[0, antenna1, antenna2, chan, 0] = self.baseline_data[bl][chan]
                        bl = bl + 1
