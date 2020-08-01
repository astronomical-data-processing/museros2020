__all__ = ['MuserBase',
           'MuserFrameHeader']

import logging
import codecs, string, struct, string
from astropy.time import Time

log = logging.getLogger('logger')

class MuserFrameHeader(object):
    """ Muser Frame Header Description Class

    """
    def __init__(self):
        self.system_time = None
        self.frequency_switch = 0
        self.bandwidth = 0
        self.quantization_leve = 0
        self.delay_switch = 0
        self.strip_switch = 0
        self.sub_band_switch = 0

class MuserBase(object):
    """ Muser Data Description Class


    """

    def __init__(self, sub_array=1):
        """

        :param sub_array: 1 - Muser Low, 2 - Muser High
        """
        self.debug = 0
        self.sub_array = sub_array
        self.total_polarization = 2
        if self.sub_array == 1:
            self.total_band = 4
        else:
            self.total_band = 32

        self.input_file_name = ''

    def convert_time(self, stime):

        # print '%x' % stime
        nanosecond = (stime & 0x3f)
        if nanosecond >= 50:
            nanosecond = 0
        nanosecond *= 20
        stime >>= 6
        # read microsecond, 6-15
        microsecond = (stime & 0x3ff)
        stime >>= 10
        # read millisecond 16-25
        millisecond = (stime & 0x3ff)
        stime >>= 10
        # read second, 26-31
        second = (stime & 0x3f)
        stime >>= 6
        # read minute, 32-37
        minute = (stime & 0x3f)
        stime >>= 6
        # read hour, 38-42
        hour = (stime & 0x1f)
        stime >>= 5
        # read day
        day = (stime & 0x1f)
        stime >>= 5
        # read month, 48-51
        month = (stime & 0x0f)
        stime >>= 4
        # read year
        year = (stime & 0xfff) + 2000
        # print tmp
        tmp = "{:4d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:d}".format(year,month,day,hour,minute,second,millisecond*1000+nanosecond)
        return Time(tmp,format='isot',scale='utc')

    def signed(self, a):
        if a >= 128 * 256 * 256:
            a = a - 256 * 256 * 256
        return a

    def convert_cross_correlation(self, buff):
        # read imaginary part of second channel
        c1 = c2 = complex()
        imag2 = struct.unpack('B', bytes([buff[8]]))[0]
        imag2 <<= 8
        imag2 |= struct.unpack('B', bytes([buff[4]]))[0]
        imag2 <<= 8
        imag2 |= struct.unpack('B', bytes([buff[0]]))[0]
        imag2 = self.signed(imag2)

        # read real part of second channel
        real2 = struct.unpack('B', bytes([buff[9]]))[0]
        real2 <<= 8
        real2 |= struct.unpack('B', bytes([buff[5]]))[0]
        real2 <<= 8
        real2 |= struct.unpack('B', bytes([buff[1]]))[0]
        real2 = self.signed(real2)

        c2 = complex(real2, imag2)

        # read imaginary part of second channel
        imag1 = struct.unpack('B', bytes([buff[10]]))[0]
        imag1 <<= 8
        imag1 |= struct.unpack('B', bytes([buff[6]]))[0]
        imag1 <<= 8
        imag1 |= struct.unpack('B', bytes([buff[2]]))[0]
        imag1 = self.signed(imag1)

        # read real part of second channel
        real1 = struct.unpack('B', bytes([buff[11]]))[0]
        real1 <<= 8
        real1 |= struct.unpack('B', bytes([buff[7]]))[0]
        real1 <<= 8
        real1 |= struct.unpack('B', bytes([buff[3]]))[0]
        real1 = self.signed(real1)

        c1 = complex(real1, imag1)
        return c1, c2

    def convert_auto_correlation(self, buff):
        r1 = struct.unpack('B', bytes([buff[3]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[2]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[1]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[0]]))[0]
        # r1 = self.signed(r1)

        r2 = struct.unpack('B', bytes([buff[7]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[6]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[5]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[4]]))[0]
        # r2 = self.signed(r2)

        return r1, r2

    def convert_time_offset(self, buff):
        r = []
        r1 = struct.unpack('B', bytes([buff[1]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[0]]))[0]

        r2 = struct.unpack('B', bytes([buff[9]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[8]]))[0]

        r.append((r2 + r1 * 0.0001) * 1e-09)

        r1 = struct.unpack('B', bytes([buff[3]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[2]]))[0]

        r2 = struct.unpack('B', bytes([buff[11]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[10]]))[0]

        r.append((r2 + r1 * 0.0001) * 1e-09)

        r1 = struct.unpack('B', bytes([buff[5]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[4]]))[0]

        r2 = struct.unpack('B', bytes([buff[13]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[12]]))[0]

        r.append((r2 + r1 * 0.0001) * 1e-09)

        r1 = struct.unpack('B', bytes([buff[7]]))[0]
        r1 <<= 8
        r1 |= struct.unpack('B', bytes([buff[6]]))[0]

        r2 = struct.unpack('B', bytes([buff[15]]))[0]
        r2 <<= 8
        r2 |= struct.unpack('B', bytes([buff[14]]))[0]

        r.append((r2 + r1 * 0.0001) * 1e-09)

        return r

    def convert_time_offset_high(self, buff):
        r = []
        r0 = struct.unpack('B', bytes([buff[0]]))[0]
        r1 = struct.unpack('B', bytes([buff[1]]))[0]
        r2 = struct.unpack('B', bytes([buff[2]]))[0]
        r3 = struct.unpack('B', bytes([buff[3]]))[0]
        r4 = struct.unpack('B', bytes([buff[4]]))[0]
        r5 = struct.unpack('B', bytes([buff[5]]))[0]
        r6 = struct.unpack('B', bytes([buff[6]]))[0]
        r7 = struct.unpack('B', bytes([buff[7]]))[0]
        r8 = struct.unpack('B', bytes([buff[8]]))[0]
        r9 = struct.unpack('B', bytes([buff[9]]))[0]
        r10 = struct.unpack('B', bytes([buff[10]]))[0]
        r11 = struct.unpack('B', bytes([buff[11]]))[0]
        r12 = struct.unpack('B', bytes([buff[12]]))[0]
        r13 = struct.unpack('B', bytes([buff[13]]))[0]
        r14 = struct.unpack('B', bytes([buff[14]]))[0]
        r15 = struct.unpack('B', bytes([buff[15]]))[0]  # read 16*8bits: 4*15(integer part) + 4*17(decimal part)

        a1 = (r15 << 7) | (r14 >> 1)
        a2 = ((r14 & 1) << 14) | (r13 << 6) | ((r12 & 0xfc) >> 2)
        a3 = ((r12 & 0x03) << 13) | (r11 << 5) | ((r10 & 0xf8) >> 3)
        a4 = ((r10 & 0x07) << 12) | (r9 << 4) | ((r8 & 0xf0) >> 4)

        b1 = ((r8 & 0x0f) << 13) | (r7 << 5) | ((r6 & 0xf1) >> 3)
        b2 = ((r6 & 0x07) << 14) | (r5 << 6) | ((r4 & 0xfc) >> 2)
        b3 = ((r4 & 0x03) << 15) | (r3 << 7) | ((r2 & 0xfe) >> 1)
        b4 = ((r2 & 0x01) << 16) | (r1 << 8) | r0

        # print b1, b2, b3, b4

        r.append((a1 + b1 * 0.00001))  # ns
        r.append((a2 + b2 * 0.00001))  # * 1e-09)
        r.append((a3 + b3 * 0.00001))  # * 1e-09)
        r.append((a4 + b4 * 0.00001))  # * 1e-09)
        return r
