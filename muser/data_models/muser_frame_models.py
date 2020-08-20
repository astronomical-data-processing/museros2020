""" MUSER-low/high Frame  Model


"""

__all__ = ['MuserFrame']

import logging
import numpy
import math
import struct

from muser.data_models.muser_data_models import MuserBase, MuserFrameHeader
from muser.data_models.muser_delay_models import MuserDelay
from muser.data_models.parameters import IF_BANDWIDTH, LOOP_MODE_LOW, NON_LOOP_MODE_LOW, NON_LOOP_MODE_HIGH, \
    LOOP_MODE_HIGH
from astropy.time import Time, TimezoneInfo
import astropy.units as u
from muser.data_models.parameters import muser_path, muser_data_path

log = logging.getLogger('logger')


class MuserFrame(MuserBase):
    def __init__(self, sub_array=1, mode=True):
        '''Construction Function

        :param sub_array 1-Muser1, 2-Muser2
        :return
        '''

        super(MuserFrame, self).__init__()
        # public member
        self.real_sub_array = 0
        self.frame_number = 8
        self.sub_channels = 0
        self.channels = 0
        self.antennas = 0
        self.dr_output_antennas = 0
        self.polarization = 0
        self.freqid = 1

        self.groups = 0
        self.fre_rows = 0
        self.source_rows = 0

        self.if_bandwidth = 25000000

        self.baseline = []

        self.channel_group = 0
        self.current_frame_time = 0
        self.current_frame_utc_time = 0
        self.polarization_number = 2
        self.first_frame_time = 0

        self.frequency = 1600000000
        self.current_frame_header = MuserFrameHeader()
        self.is_loop_mode = mode
        self.read_number = 0

        # ntimes x antenna x antenna x channels x polarization
        self.is_data_buffer_initialized = False
        # Define sub array: 1 - Muser Low, 2 - Muser High
        self.set_array(sub_array)

    def set_array(self, sub_array=1):
        '''
        switch between muser-i and muser-ii
        sub_array: 1:muser-i/2:muser-ii
        '''

        muser_para = {0: (16, 40, 44, 4, 2), 1: (16, 60, 64, 33, 2)}

        self.sub_array = sub_array
        if self.sub_array not in [1, 2]:
            return False

        self.sub_channels, self.antennas, self.dr_output_antennas, self.frame_number, self.polarization_number = \
            muser_para[self.sub_array - 1]

        log.debug('Array: muser-%d: antennas:%d, frame: %d, polarization:%d' % (self.sub_array,
                                                                                self.antennas, self.frame_number,
                                                                                self.polarization_number))

        self.correlation_data = numpy.zeros(
            shape=(self.sub_channels, self.dr_output_antennas, self.dr_output_antennas),
            dtype=float)

        self.par_delay = numpy.ndarray(shape=(self.dr_output_antennas), dtype=float)

        self.delay = numpy.ndarray(shape=(self.dr_output_antennas), dtype=float)

        self.delay_compensation = MuserDelay()

        # self.env.antenna_loaded = False

    def init_data_buffer(self):
        self.block_data = numpy.zeros((self.antennas, self.antennas, self.sub_channels), dtype=complex)
        if self.is_loop_mode:
            self.block_full_data = numpy.zeros(
                (self.antennas, self.antennas, self.sub_channels * self.frame_number, self.polarization_number),
                dtype=complex)
            self.real_polarization_number = 2
            self.real_frame_number = self.frame_number
            if self.sub_array == 1:
                self.start_frequency = 400
                self.end_frequency = 2000
            else:
                self.start_frequency = 2000
                self.end_frequency = 15000
        else:
            self.block_full_data = numpy.zeros((self.antennas, self.antennas, self.sub_channels, 1), dtype=complex)
            self.real_polarization_number = 1
            self.real_frame_number = 1
            if self.sub_array ==1:
                self.start_frequency = self.sub_band * 400
                self.end_frequency = self.sub_band * 400 + 400
            else:
                if self.sub_band == 32:
                    self.start_frequency = 14600
                    self.end_frequency = 15000
                self.start_frequency = self.sub_band * 400 + 2000
                self.end_frequency = (self.sub_band + 1) * 400 + 2000

        self.is_data_buffer_initialized = True

    def search_frame_header_with_byte(self):
        # search the header of each frame which sent by digital receiver
        # per frame/3ms
        # in_file: the file handle of the raw file
        # return True : find a valid frame

        while True:
            bb = self.in_file.read(1)
            if not bb:
                break
            beginner = struct.unpack('b', bb)[0]
            # print '%2x ' % beginner

            if beginner == 0x55:
                for i in range(0, 27):
                    bb = self.in_file.read(1)
                    if not bb: break
                    second = struct.unpack('b', bb)[0]
                    # print ('%d %2x ') % (i,second)
                    if second != 0x55:
                        break
                if i == 26:
                    self.in_file.seek(4, 1)
                    return True
        return False

    def search_frame_header(self):
        # search the header of each frame which sent by digital receiver
        # per frame/3ms
        # in_file: the file handle of the raw file
        # return True : find a valid frame
        offset = [100000, 204800]
        pos = self.in_file.tell()
        if pos > 0:
            if pos % offset[self.sub_array - 1] != 0:
                self.in_file.seek(((pos // offset[self.sub_array - 1]) + 1) * offset[self.sub_array - 1])
        self.in_file.seek(32, 1)
        if self.debug:
            log.debug('Frame header searched - %d ', (pos))

        return True

    def seek_frame_header(self):
        # search the header of each frame which sent by digital receiver
        # per frame/3ms
        # in_file: the file handle of the raw file
        # return True : find a valid frame
        offset = [100000, 204800]
        pos = self.in_file.tell()
        if pos > 0:
            self.in_file.seek((pos // offset[self.sub_array - 1]) * offset[self.sub_array - 1])
        self.in_file.seek(32, 1)
        return True

    def search_next_full_frame_head(self):
        while True:
            if self.read_one_frame() == False:
                self.in_file.close()
                return False
            if self.is_loop_mode == False:
                break
            else:
                if self.is_loop_mode and self.sub_band == 0 and self.polarization == 0:  # Find file in previous 1 minute
                    break
        self.skip_frames(0)
        log.debug('Frame located.')
        return True

    def skip_frames(self, number_of_frames):
        offset = [100000, 204800]
        pos = self.in_file.tell()
        if pos > 0 and (pos % offset[self.sub_array - 1] != 0):
            self.in_file.seek((pos // offset[self.sub_array - 1]) * offset[self.sub_array - 1])

        self.in_file.seek(offset[self.sub_array - 1] * number_of_frames, 1)
        log.debug('Skip %d frames' % (number_of_frames))
        return True

    def read_one_frame(self, search=True):
        if search == True:
            if (self.search_frame_header() == False):
                log.error('Cannot find a correct header information.')
                return False
        else:
            self.in_file.seek(32, 1)
        tmp = self.in_file.read(8)
        tmp_time = struct.unpack('Q', tmp)[0]

        self.current_frame_time = self.convert_time(tmp_time)
        # if self.Debug:
        log.debug("Read frame date and time: %s " % self.current_frame_time.isot)

        self.current_frame_date_time = self.current_frame_time

        # the date and time are beijing time of china, utc = cst - 8
        # utc = cst - 8
        self.current_frame_utc_time = self.current_frame_time - 8 * u.hour

        if (self.current_frame_time < Time('2015-01-01T00:00:00', format='isot')):
            self.version = False
        else:
            self.version = True

        # read out parameters
        self.current_frame_header.frequency_switch = struct.unpack('I', self.in_file.read(4))[0]
        self.in_file.seek(128, 1)
        self.current_frame_header.bandwidth = struct.unpack('I', self.in_file.read(4))[0]
        self.current_frame_header.quantization_level = struct.unpack('I', self.in_file.read(4))[0]
        self.current_frame_header.delay_switch = struct.unpack('I', self.in_file.read(4))[0]
        self.current_frame_header.strip_switch = struct.unpack('I', self.in_file.read(4))[0]
        self.current_frame_header.sub_band_switch = struct.unpack('I', self.in_file.read(4))[0]
        self.if_bandwidth = IF_BANDWIDTH[self.current_frame_header.bandwidth]

        if self.version == True:
            if (self.current_frame_header.sub_band_switch & 0xffff) == 0x3333:  # Polarization
                self.polarization = 0
            else:
                self.polarization = 1
        else:
            self.polarization = (self.read_number % 2)
            self.read_number = self.read_number + 1
        self.sub_band_switch = self.current_frame_header.sub_band_switch >> 16  # Sub_band frequency
        self.is_loop_mode = False

        if self.sub_array == 1:
            # We have to support two types of files
            if self.version == False:
                if self.sub_band_switch in [0x3333, 0x7777, 0xbbbb, 0xcccc]:
                    self.channel_group = NON_LOOP_MODE_LOW[self.sub_band_switch][0]
                    self.frequency = NON_LOOP_MODE_LOW[self.sub_band_switch][1]
                    self.freqid = NON_LOOP_MODE_LOW[self.sub_band_switch][2]

                    self.numrows = 1
                    self.groups = 60 * 1000 / 25 * 8 * (self.antennas * (self.antennas - 1) // 2)
                    self.source_rows = 60 * 1000 / 25 * 8 / 2

                else:
                    self.is_loop_mode = True

                    self.numrows = 4
                    self.groups = 60 * 1000 / 25 * 8 * (self.antennas * (self.antennas - 1) // 2) / 2
                    self.source_rows = 60 * 1000 / 25 * 8

                    # In loop mode, we should read another two bytes in absolute offset: 99264
                    # We have read 192 bytes. So, we need to further move forward 99264-192  bytes
                    # 0000,5555,AAAA,FFFF 0.4-0.8GHZ 0.8~1.2GHz 1.2~1.6GHz 1.6~2.0GHz
                    self.in_file.seek(99264 - 192, 1)
                    self.sub_band_switch = struct.unpack('H', self.in_file.read(2))[0]
                    if self.sub_band_switch == 0x0000:
                        self.channel_group = 0
                        self.frequency = 400000000
                        self.freqid = 1
                    elif self.sub_band_switch == 0x5555:
                        self.channel_group = 16
                        self.frequency = 800000000
                        self.freqid = 2
                    elif self.sub_band_switch == 0xaaaa:
                        self.channel_group = 32
                        self.frequency = 1200000000
                        self.freqid = 3
                    elif self.sub_band_switch == 0xffff:
                        self.channel_group = 48
                        self.frequency = 1600000000
                        self.freqid = 4
                    self.in_file.seek(-(99264 - 192 + 2), 1)
            else:
                if self.sub_band_switch in NON_LOOP_MODE_LOW:
                    self.channel_group = NON_LOOP_MODE_LOW[self.sub_band_switch][0]
                    self.frequency = NON_LOOP_MODE_LOW[self.sub_band_switch][1]
                    self.freqid = NON_LOOP_MODE_LOW[self.sub_band_switch][2]
                else:
                    self.is_loop_mode = True
                    # in loop mode, we should read another two bytes in absolute offset: 99264
                    # we have read 192 bytes. so, we shoudl foward 99264-192  bytes
                    # 0000,5555,aaaa,ffff 0.4-0.8ghz 0.8~1.2ghz 1.2~1.6ghz 1.6~2.0ghz
                    self.in_file.seek(16, 1)  # todo:
                    self.current_frame_header.sub_band_switch = struct.unpack('B', self.in_file.read(1))[0]
                    # print "self.current_frame_header.sub_band_switch", self.current_frame_header.sub_band_switch
                    self.polarization = LOOP_MODE_LOW[self.current_frame_header.sub_band_switch][0]
                    self.channel_group = LOOP_MODE_LOW[self.current_frame_header.sub_band_switch][1]
                    self.frequency = LOOP_MODE_LOW[self.current_frame_header.sub_band_switch][2]
                    self.freqid = LOOP_MODE_LOW[self.current_frame_header.sub_band_switch][3]
                    self.in_file.seek(-17, 1)
                self.sub_band = self.channel_group // 16

        elif self.sub_array == 2:
            if self.sub_band_switch in NON_LOOP_MODE_HIGH:  # Sub_band frequency
                if self.sub_band_switch not in NON_LOOP_MODE_HIGH:
                    self.sub_band_switch = 0
                self.channel_group = NON_LOOP_MODE_HIGH[self.sub_band_switch][0]
                self.frequency = NON_LOOP_MODE_HIGH[self.sub_band_switch][1]
                self.freqid = NON_LOOP_MODE_HIGH[self.sub_band_switch][2]
            else:
                self.is_loop_mode = True
                self.in_file.seek(4, 1)
                self.current_frame_header.sub_band_switch = struct.unpack('H', self.in_file.read(2))[0]
                self.polarization = LOOP_MODE_HIGH[self.current_frame_header.sub_band_switch][0]
                self.channel_group = LOOP_MODE_HIGH[self.current_frame_header.sub_band_switch][1]
                self.frequency = LOOP_MODE_HIGH[self.current_frame_header.sub_band_switch][2]
                self.freqid = 1
                self.in_file.seek(-6, 1)

        self.sub_band = self.channel_group // 16
        if self.is_loop_mode and self.sub_band == 0 and self.polarization == 0:
            self.first_frame_time = self.current_frame_time
        elif self.is_loop_mode == False:
            self.first_frame_time = self.current_frame_time

        if self.version == True:
            obs_target = struct.unpack('H', self.in_file.read(2))[0]  # calibration source
            if obs_target == 0x0000:
                self.obs_target = "sun"
            elif obs_target == 0x0101:
                self.obs_target = "satellite"
            elif obs_target == 0x0202:
                self.obs_target = "swan"
            else:
                self.obs_target = "other"

            tag = struct.unpack('H', self.in_file.read(2))[0]

            if tag == 0x0000:  # high array or low array
                self.real_sub_array = 1
            elif tag == 0x4e4e:
                self.real_sub_array = 2
            if self.sub_array == 1:
                self.in_file.seek(-4, 1)
        else:
            self.obs_target = "sun"
            self.real_sub_array = self.sub_array

        # To simply programming, we use one variable to support Numpy Array operation
        if self.is_loop_mode:
            self.real_polarization = self.polarization
            self.real_sub_band = self.sub_band
        else:
            self.real_polarization = 0
            self.real_sub_band = 0
        # Initialized data reading buffer for first time ONCE
        if not self.is_data_buffer_initialized:
            self.init_data_buffer()
        log.debug(
            'Read header info - time: %s - array:%d band:%d polization:%d frequence:%d' % (
                self.current_frame_time.isot, self.real_sub_array, self.sub_band, self.polarization,
                self.frequency))
        return True

    def read_data(self):
        if self.sub_array == 1:
            log.debug('Read data current pos %d' % self.in_file.tell())
            self.in_file.seek(2752, 1)

            for channel in range(0, self.sub_channels, 2):
                bl = 0
                for antenna1 in range(0, self.dr_output_antennas - 1):
                    for antenna2 in range(antenna1 + 1, self.dr_output_antennas):
                        buff = self.in_file.read(12)
                        c1, c2 = self.convert_cross_correlation(buff)
                        if (antenna1 < self.antennas and antenna2 < self.antennas ):
                            self.block_data[antenna2, antenna1, channel] = c1
                            self.block_data[antenna2, antenna1, channel + 1] = c2
                            bl = bl + 1
                self.in_file.seek(40, 1)

                for antenna in range(0, self.dr_output_antennas, 4):

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    if antenna < self.antennas - 4:
                        self.block_data[antenna, antenna, channel] = r1
                        self.block_data[antenna, antenna, channel + 1] = r2

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    if antenna < self.antennas - 4:
                        self.block_data[antenna + 1, antenna + 1, channel] = r1
                        self.block_data[antenna + 1, antenna + 1, channel + 1] = r2

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    if antenna < self.antennas - 4:
                        self.block_data[antenna + 2, antenna + 2, channel] = r1
                        self.block_data[antenna + 2, antenna + 2, channel + 1] = r2

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    if antenna < self.antennas - 4:
                        self.block_data[antenna + 3, antenna + 3, channel] = r1
                        self.block_data[antenna + 3, antenna + 3, channel + 1] = r2
                    self.in_file.seek(32, 1)

                    if channel == self.sub_channels - 2:
                        self.in_file.seek(-24, 1)
                        buff = self.in_file.read(16)
                        (self.par_delay[antenna], self.par_delay[antenna + 1], self.par_delay[antenna + 2],
                         self.par_delay[antenna + 3]) = self.convert_time_offset(buff)[:]
                        self.in_file.seek(8, 1)

                # print self.auto_correlation_data
                self.in_file.seek(32, 1)

        elif self.sub_array == 2:
            self.in_file.seek(1692, 1)  # 1690 reserve and 2bytes frequency code
            visbility = numpy.zeros(shape=(self.dr_output_antennas * (self.dr_output_antennas - 1) / 2),
                                    dtype=complex)
            for channel in range(0, self.sub_channels):  # read data of all antennas in one channel
                for bl_len in range(0, self.dr_output_antennas * (self.dr_output_antennas - 1) / 2, 2):
                    buff = self.in_file.read(12)
                    c1, c2 = self.convert_cross_correlation(buff)
                    visbility[bl_len] = c2
                    visbility[bl_len + 1] = c1

                bl1, bl2 = 0, 0
                for antenna1 in range(0, self.antennas - 1):
                    for antenna2 in range(antenna1 + 1, self.dr_output_antennas):
                        if antenna2 < self.antennas:
                            self.block_data[antenna2, antenna1, channel] = visbility[bl2]
                            bl1 += 1
                        bl2 += 1

                newch = channel % 2
                for antenna in range(0, 8):
                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    self.block_data[0 + antenna * 8 + newch * 4, 0 + antenna * 8 + newch * 4, 0 + channel - newch] = r1
                    self.block_data[0 + antenna * 8 + newch * 4, 0 + antenna * 8 + newch * 4, 1 + channel - newch] = r2

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    self.block_data[1 + antenna * 8 + newch * 4, 1 + antenna * 8 + newch * 4, 0 + channel - newch] = r1
                    self.block_data[1 + antenna * 8 + newch * 4, 1 + antenna * 8 + newch * 4, 1 + channel - newch] = r2

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    self.block_data[2 + antenna * 8 + newch * 4, 2 + antenna * 8 + newch * 4, 0 + channel - newch] = r1
                    self.block_data[2 + antenna * 8 + newch * 4, 2 + antenna * 8 + newch * 4, 1 + channel - newch] = r2

                    buff1 = self.in_file.read(8)
                    r1, r2 = self.convert_auto_correlation(buff1)
                    self.block_data[3 + antenna * 8 + newch * 4, 3 + antenna * 8 + newch * 4, 0 + channel - newch] = r1
                    self.block_data[3 + antenna * 8 + newch * 4, 3 + antenna * 8 + newch * 4, 1 + channel - newch] = r2

                    self.in_file.seek(32, 1)
                    if channel == 14 or channel == 15:
                        self.in_file.seek(-24, 1)
                        buff = self.in_file.read(16)
                        (self.par_delay[0 + antenna * 8 + newch * 4], self.par_delay[1 + antenna * 8 + newch * 4],
                         self.par_delay[2 + antenna * 8 + newch * 4],
                         self.par_delay[3 + antenna * 8 + newch * 4]) = self.convert_time_offset_high(buff)[:]

                        self.in_file.seek(8, 1)
                self.in_file.seek(72, 1)

        log.debug("Read visibility  and auto correlation data.")
        return True

    def read_one_file(self):  # data in the same channe: 780*2400
        if self.search_first_frame() == False:
            log.error("Cannot find observational data.")
            return str(False), []
        self.seek_frame_header()
