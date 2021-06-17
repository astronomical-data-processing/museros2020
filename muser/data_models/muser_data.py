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
from muser.components.ephem.sun_position import get_sun

log = logging.getLogger('muser')


class MuserData(MuserFrame):
    def __init__(self, sub_array=1, mode=True, file_name=None, start_time=None):

        super(MuserData, self).__init__(sub_array=sub_array, mode=mode)
        if file_name is not None:
            self.input_file_name = file_name
        self.if_time_setup = False
        self.if_file_opened = False
        self.in_file = None
        self.integral_time = 0
        self.init_data = False

        self.last_calibration_priority = -1
        self.last_date_time = 0
        self.last_sub_array = -1

        self.current_file_name = ''
        self.file_list = []
        self.if_file_opened = False
        self.need_change_file = False
        if file_name is not None:
            self.input_file_name = file_name
        if start_time is not None:
            self.start_date_time = Time(start_time, format='isot')
            # self.first_frame_time = Time(start_time, format='isot')

    def __del__(self):
        self.close_file()

    # @property
    # def filename(self, name):
    #     self.filename = name

    def set_date_time(self, time):
        self.start_date_time = Time(time, format='isot')
        self.first_frame_time = Time(time, format='isot')

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

    def search_first_file(self, frame_time=None):
        '''
        search the first file of processing
        '''
        if len(self.input_file_name) == 0:
            if frame_time is None:
                return False
            self.start_date_time = Time(frame_time, format='isot')

        if self.open_data_file() == False:
            log.debug("Cannot open observational data.")
            return False

        # Read first frame and check the observational date and tim  e
        if self.read_one_frame() == False:
            log.debug("Cannot read a frame.")
            return False

        # if self.current_file_name != '':  # Start!!!
        #     self.start_date_time = self.current_frame_time

        if self.start_date_time < self.current_frame_time:
            # previous 1 minute
            log.debug("Change observational file.")
            if self.open_next_file(-1) == False:
                return False
            # Read first frame and check the observational date and time
            if self.read_one_frame() == False:
                self.in_file.close()
                return False
        # self.skip_frames(0)
        log.debug("Found observational data file: %s" % (os.path.basename(self.current_file_name)))
        return True

    def open_next_file(self, time_minute=1):
        search_date_time = self.current_frame_time + time_minute * u.minute

        full_file_name = self.muser_data_file_name(search_date_time.datetime.year,
                                                   search_date_time.datetime.month,
                                                   search_date_time.datetime.day,
                                                   search_date_time.datetime.hour,
                                                   search_date_time.datetime.minute)
        self.input_file_name = full_file_name
        return self.open_data_file()

    def muser_data_file_name(self, year, month, day, hour, minute):
        file_name = ('CSRH_%04d%02d%02d-%02d%02d') % (year, month, day, hour, minute)
        for file in self.file_list:
            if file_name in file:
                return muser_data_path(file)
        return ''
        # print("Cannot find observational data or not a MUSER file.")
        # exit(1)

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

    def close_file(self):
        try:
            if self.if_file_opened:
                self.in_file.close()
        finally:
            pass

    def open_data_file(self):
        if_input_file_name = False
        if self.input_file_name == '':
            full_file_name = self.muser_data_file_name(self.start_date_time.datetime.year,
                                                       self.start_date_time.datetime.month,
                                                       self.start_date_time.datetime.day,
                                                       self.start_date_time.datetime.hour,
                                                       self.start_date_time.datetime.minute)
        else:
            full_file_name = self.input_file_name
            if_input_file_name = True

        if not os.path.exists(full_file_name):
            log.error("Cannot find file: %s " % (full_file_name))
            return False
        try:
            self.in_file = open(full_file_name, 'rb')
            self.current_file_name = full_file_name
            self.if_file_opened = True
            offset = [100000, 204800]
            # Get last frame
            self.in_file.seek(0, 0)
            self.in_file.seek(-offset[self.sub_array - 1], 2)
            self.read_one_frame()
            self.file_end_time = self.current_frame_time
            # Get first frame
            self.in_file.seek(0, 0)
            self.read_one_frame()
            self.file_first_time = self.current_frame_time
            if if_input_file_name:
                self.start_date_time = self.current_frame_time
            #  Reset file pointer
            self.in_file.seek(0, 0)
            log.debug("File opened: %s" % (os.path.basename(self.current_file_name)))
            log.debug('First and Last frame time: %s - %s' % (self.file_first_time.isot, self.file_end_time.isot))

            return True
        except IOError:
            # self.in_file.close()
            return False

    def search_frame(self, search_time):
        '''
        search first frame with proper date and time
        '''
        self.in_file.seek(0, 0)
        if not self.read_one_frame():
            return False
        self.start_date_time = Time(search_time, format='isot')
        if search_time > self.file_end_time:
            if not self.open_next_file(1):
                return False
        if search_time < self.file_first_time:
            if not self.open_next_file(-1):
                return False

        if search_time < self.current_frame_time:
            t_offset = self.current_frame_time - self.start_date_time
        else:
            t_offset = self.start_date_time - self.current_frame_time

        log.debug('Current frame time: %s' % (self.current_frame_time.isot))

        # Calculate time offset with unit of microsecond
        time_offset = t_offset.to_value(unit='microsecond')

        # time_offset = t_offset.datetime.second * 1e6 + t_offset.datetime.microsecond

        skip_frame_number = int(time_offset / 3125)
        if skip_frame_number >= 20:
            skip_frame_number = skip_frame_number - 20

        log.debug('Time interval %d, skip frames: %d' % (time_offset, skip_frame_number))

        if (skip_frame_number >= 2):
            self.skip_frames(skip_frame_number)

        # If can find, search first frame
        while True:
            if self.read_one_frame() == False:
                return False
            if self.is_loop_mode == False:
                if self.start_date_time <= self.current_frame_time:
                    break
            else:
                if search_time >= self.current_frame_time and self.sub_band == 0 and self.polarization == 0:  # Find file in previous 1 minute
                    break
            if self.current_frame_time == self.file_end_time:
                if not self.open_next_file(1):
                    return False

        log.debug('Frame located.')
        self.skip_frames(0)
        return True

    def search_first_frame(self):
        '''
        search first frame with proper date and time
        '''

        if self.search_first_file(frame_time=self.start_date_time) == False:
            log.error("Search first frame: Cannot find a proper frame from the observational data.")
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
                if self.is_loop_mode and self.start_date_time <= self.current_frame_time \
                        and self.sub_band == 0 and self.polarization == 0:  # Find file in previous 1 minute
                    self.first_frame_time = self.current_frame_time
                    break
                if self.is_loop_mode == False and self.start_date_time <= self.current_frame_time:
                    self.first_frame_time = self.current_frame_time
                    break

            if self.read_one_frame() == False:
                self.in_file.close()
                return False
        log.debug('Frame located.')
        return True

    def read_one_data_for_full(self, ):
        self.read_data()
        from copy import deepcopy
        self.block_full_data[:, :,
        self.real_sub_band * self.sub_channels: self.real_sub_band * self.sub_channels + 16,
        self.real_polarization] = deepcopy(self.block_data[:, :, :])

    def read_full_frame(self, search=True, read_data=False):
        index = 1
        while True:
            self.block_full_data *= 0.
            total_frames = self.frame_number * self.polarization_number
            if self.is_loop_mode:
                frame = 0
                current_time = self.current_frame_time
                while frame < total_frames:
                    log.info("Reading No. %d %s %d %d" % (
                        frame, self.current_frame_time.isot, self.sub_band,
                        self.polarization))
                    if not self.read_one_frame():
                        return False
                    log.debug("Frame time:{} {} {}".format(self.current_frame_time, self.sub_band,self.polarization))
                    if read_data:
                        self.read_one_data_for_full()
                    # print(self.current_frame_time, current_time, self.sub_band, self.polarization,
                    #       (self.current_frame_time - self.first_frame_time).to_value('s'))
                    if (self.current_frame_time - current_time).to_value('s') >= 4 / 1000.:
                        log.info("Find frame lost or missing.")
                        if self.search_frame(self.current_frame_time.isot):
                            break
                        else:
                            return False
                    else:
                        frame = frame + 1
                        current_time = self.current_frame_time
                    if self.current_frame_time == self.file_end_time:
                        if not self.open_next_file(1):
                            return False
                else:
                    break
            else:
                if self.current_frame_time == self.file_end_time:
                    if not self.open_next_file(1):
                        return False
                if not self.read_one_frame():
                    return False
                if read_data:
                    self.read_one_frame_for_full()
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

    # def delay_process(self, planet):
    #     # print "delay processing..."
    #     parameter = 0.
    #     delay = numpy.zeros(shape=(self.dr_output_antennas), dtype=float)
    #
    #     if self.sub_array == 1:  # muser-1
    #         if planet == 'sun':
    #             parameter = 12.5
    #         elif planet == 'satellite':
    #             parameter = 2.5
    #         delayns = self.delay_compensation.get_delay_value(self.sub_array,
    #                                                           self.current_frame_time)
    #         delay = self.par_delay * (10 ** 9) - delayns
    #     else:  # muser-2
    #         parameter = 12.5
    #         delay = self.par_delay
    #
    #     delay_real = delay[0:39]
    #     [delay_x, delay_y] = numpy.meshgrid(delay_real, delay_real)
    #     if self.sub_array == 1:
    #         delay_matrix = delay_x - delay_y
    #         delay_matrix_int = numpy.trunc(delay_x) - numpy.trunc(delay_y)
    #
    #         frequency = (numpy.arange(400, 2000, 25) + parameter) / 1000.
    #         frequency_interval = (numpy.arange(0, self.sub_channels) * 25 + parameter + 50) / 1000.
    #
    #     for channel in range(0, self.sub_channels):
    #         bl = 0
    #         for antenna1 in range(0, self.antennas - 1):  # SubChannelsLow = 16
    #             for antenna2 in range(antenna1 + 1, self.antennas):
    #                 tg = delay[antenna2] - delay[antenna1]
    #                 tg0 = int(delay[antenna2]) - int(delay[antenna1])
    #                 if self.sub_array == 1:
    #                     Frf = (self.frequency * 1e-6 + channel * 25 + parameter) / 1000.0
    #                     Fif = (channel * 25 + parameter + 50.0) / 1000.0
    #                     phai = 2 * numpy.pi * (Frf * tg - Fif * tg0)
    #                     self.block_data[antenna1, antenna2, channel] = complex(
    #                         self.block_data[antenna1, antenna2, channel].real * numpy.cos(phai) +
    #                         self.block_data[antenna1, antenna2, channel].imag * numpy.sin(phai),
    #                         self.block_data[antenna1, antenna2, channel].imag * numpy.cos(phai) -
    #                         self.block_data[antenna1, antenna2, channel].real * numpy.sin(phai))
    #                     # self.baseline_data[bl][channel] = complex(
    #                     #     self.baseline_data[bl][channel].real * numpy.cos(phai) +
    #                     #     self.baseline_data[bl][channel].imag * numpy.sin(phai),
    #                     #     self.baseline_data[bl][channel].imag * numpy.cos(phai) -
    #                     #     self.baseline_data[bl][channel].real * numpy.sin(phai))
    #                 else:
    #                     Frf = (self.frequency * 1e-6 + (15 - channel) * 25 + parameter) / 1000.0
    #                     Fif = (channel * 25 + parameter + 50.0) / 1000.0  # local frequency(GHz)
    #                     phai = 2 * numpy.pi * (-Frf * tg - Fif * tg0)
    #                     # phai = 2 * pi * Fif * tg0 + 2 * pi * Frf * (tg - tg0)
    #                     self.block_data[antenna1, antenna2, channel] = complex(
    #                         self.block_data[antenna1, antenna2, channel].real * numpy.cos(phai) +
    #                         self.block_data[antenna1, antenna2, channel].imag * numpy.sin(phai),
    #                         self.block_data[antenna1, antenna2, channel].imag * numpy.cos(phai) -
    #                         self.block_data[antenna1, antenna2, channel].real * numpy.sin(phai))
    #
    #                     # self.baseline_data[bl][channel] = complex(
    #                     #     self.baseline_data[bl][channel].real * numpy.cos(phai) +
    #                     #     self.baseline_data[bl][channel].imag * numpy.sin(phai),
    #                     #     self.baseline_data[bl][channel].imag * (-1) * numpy.cos(phai) +
    #                     #     self.baseline_data[bl][channel].real * numpy.sin(phai))
    #                 bl = bl + 1
    #
    #     log.debug("Delay Process and fringe stopping... Done.")

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

        delay_real = delay[0:self.antennas]
        [delay_x, delay_y] = numpy.meshgrid(delay_real, delay_real)
        delay_matrix = delay_x - delay_y
        delay_matrix_int = numpy.ceil(delay_x) - numpy.ceil(delay_y)

        freq = (numpy.arange(self.start_frequency, self.end_frequency, 25) + parameter) / 1000.
        freq_interval = numpy.repeat((numpy.arange(0, self.sub_channels) * 25 + parameter + 50) / 1000., 4)
        phai1 = numpy.einsum('ij,kl->ijk', delay_matrix, freq.reshape(-1, 1))
        phai2 = numpy.einsum('ij,kl->ijk', delay_matrix_int, freq_interval.reshape(-1, 1))
        phai_block = 2 * numpy.pi * (phai1 - phai2)

        for pol in range(self.real_polarization_number):
            real = self.block_full_data[:, :, :, pol].real * numpy.cos(phai_block) + self.block_full_data[:, :, :,
                                                                                     pol].imag * numpy.sin(phai_block)
            imag = self.block_full_data[:, :, :, pol].imag * numpy.cos(phai_block) - self.block_full_data[:, :, :,
                                                                                     pol].real * numpy.sin(phai_block)
            self.block_full_data[:, :, :, pol] = numpy.vectorize(complex)(real, imag)
        #
        # print(phai_block[11,10,0])
        # for count in range(self.real_frame_number):
        #     for channel in range(0, self.sub_channels):
        #         bl = 0
        #         for antenna1 in range(0, self.antennas - 1):  # SubChannelsLow = 16
        #             for antenna2 in range(antenna1 + 1, self.antennas):
        #                 tg = delay[antenna2] - delay[antenna1]
        #                 tg0 = int(delay[antenna2]) - int(delay[antenna1])
        #                 if self.sub_array == 1:
        #                     if self.is_loop_mode:
        #                         Frf = (400+(count%2)*400 + channel * 25 + parameter) / 1000.0
        #                     else:
        #                         Frf = (self.frequency*1e06 + channel * 25 + parameter) / 1000.0
        #                     Fif = (channel * 25 + parameter + 50.0) / 1000.0
        #                     phai = 2 * numpy.pi * (Frf * tg - Fif * tg0)
        #                 else:
        #                     Frf = (self.frequency * 1e-6 + (15 - channel) * 25 + parameter) / 1000.0
        #                     Fif = (channel * 25 + parameter + 50.0) / 1000.0  # local frequency(GHz)
        #                     phai = 2 * numpy.pi * (-Frf * tg - Fif * tg0)
        #                 if antenna2 == 11 and antenna1 ==10 and count==0 and channel == 0:
        #                     print("Phai:",phai)
        #                 for pol in range(self.real_polarization_number):
        #                     cc = count * 16 + channel
        #                     self.block_full_data[antenna1, antenna2, count * 16 + channel, pol] = complex(
        #                         self.block_full_data[antenna2, antenna1, cc, pol].real * numpy.cos(phai) +
        #                         self.block_full_data[antenna2, antenna1, cc, pol].imag * numpy.sin(phai),
        #                         self.block_full_data[antenna2, antenna1, cc, pol].imag * numpy.cos(phai) -
        #                         self.block_full_data[antenna2, antenna1, cc, pol].real * numpy.sin(phai))
        #                 bl = bl + 1

        log.debug("Block Data Delay Process and fringe stopping... Done.")

    def count_frame_number(self, time_start, time_end):
        self.start_frame_time = Time(time_start, format='isot')
        count = 0
        end_frame_time = Time(time_end, format='isot')
        while True:
            if self.read_full_frame():
                if time_end is None:
                    break
                if self.first_frame_time > end_frame_time:
                    break
                count = count + 1
            else:
                break
        return count

    def phase_calibration(self, cal):
        log.debug("Satellite phase correction")
        # cal = cal.reshape(self.block_full_data.shape)
        if self.sub_array == 1:
            amplitude = abs(self.block_full_data)
            phai_sun = numpy.arctan2(self.block_full_data.imag, self.block_full_data.real)
            phai_sat = numpy.arctan2(cal.imag[0, ...], cal.real[0, ...])
            phai = phai_sun - phai_sat
            real = amplitude * numpy.cos(-phai)
            imag = amplitude * numpy.sin(-phai)
            self.block_full_data = numpy.vectorize(complex)(real, imag)
