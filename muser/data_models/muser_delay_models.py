"""Optical Transfer Delay Calibration


"""

import numpy
import logging

class MuserDelay(object):
    __instance = None

    def __new__(cls):
        if cls.__instance is None:
            cls.__instance = super(MuserDelay, cls).__new__(cls)
            cls.__instance.__initialized = False
        return cls.__instance

    def __init__(self):
        if (self.__initialized): return
        self.__initialized = True

        self.delay_from_text()

    def delay_from_text(self):
        self.delay_array_low = numpy.array(
            [-56, 0, 48, 921, 13, 59, -3, 460, 49, 69, -675, -157, 363, -65, 30, 42, 51, 121, -2, 73, 35,
             26, 74, 35, -3, 47, -71, 75, 343, 56, 32, 313, 678, 12, -30, 48, -18, 20, 10, -1666, 0, 0, 0,
             0], dtype="int")

        DelayShift_ns1 = numpy.array(
            [0, 0, -30, -20, 10, 9, -15, -1, -48, -33, 0, 186, 0, 0, 0, -26, 0, 0, 0, 0, -24, -1, -57, 0,
             0, 3, 66, -12, 0, -20, -1, 0, -22, 0, 7, 10, 2, 1, 0, 0, 0, 0, 0, 0], dtype="int")
        DelayShift_ns2 = numpy.array(
            [0, 0, 0, 0, 0, 0, 0, 0, 12, 10, 0, 0, 0, 0, -31, 45, 0, 0, 0, 0, 11, 20, 0, 0, 0, 0, 0, 9, 11,
             0, 10, 10, 9, 1, 10, 4, 11, 0, 0, 0, 0, 0, 0, 0], dtype="int")
        DelayShift_ns3 = numpy.array(
            [0, 0, -102, -50, -51, -47, -48, -50, -62, -64, 0, 0, 0, 0, -57, -57, -60, -60, -60, -60, -60,
             -69, -70, -54, -37, -47, 0, -60, -62, -52, -61, -59, -59, -49, -63, -65, -59, -60, -67, 0, 0,
             0, 0, 0], dtype="int")
        DelayShift_ns4 = numpy.array(
            [106, 0, 101, 50, 48, 46, 53, 50, 52, 50, 0, 0, -17, 0, 26, 54, 0, 0, 0, 0, 50, 50, 69, 46, 47,
             50, 0, 50, 50, 53, 50, 49, 49, 47, 47, 51, 48, 53, 51, 0, 0, 0, 0, 0], dtype='int')
        DelayShift_ns5 = numpy.array(
            [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
             -1, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='int')

        for i in range(0, self.delay_array_low.size):  # dr_output_antennas = 44
            self.delay_array_low[i] = self.delay_array_low[i] + DelayShift_ns1[i] + DelayShift_ns2[i] + DelayShift_ns3[
                i] + DelayShift_ns4[i] + DelayShift_ns5[i]

        self.delay_array_hi = numpy.array([0, 0, 10, -28, 0, 477, 23, -18,  # $   ;a1-a8               ;20140916
                                           0, -20, 1, 30, 0, 0, -14, 61,  # $   ;a9-a16
                                           0, -5, 40, 123, 0, -45, -15, 77,  # $   ;a17-a20,b1-b4
                                           0, 232, 40, 321, 0, -114, -17, -7,  # $   ;b5-b12
                                           0, -47, -29, 148, 0, 0, 0, 112,  # $   ;b13-b20
                                           0, -30, -43, -11, 0, 34, -14, 28,  # $   ;c1-c8
                                           0, 196, 36, 193, 0, -8, -33, -59,  # $   ;c9-c16
                                           0, -40, -2, -27, 0, 0, 0, 0], dtype="int")  # ;c17-c20,d1-d4
        delay_add1 = numpy.array([0, 0, -3, 3, -10, -3, -1, 0,  # $             ;20141023
                                  -15, 0, -2, -1, 145, 0, 3, 0,  # $
                                  -60, 0, 0, -2, -22, 0, 0, 0,  # $
                                  0, 58, 0, 0, 150, -2, -2, -1,  # $
                                  -16, -1, -3, 0, -14, 0, 0, 0,  # $
                                  0, 0, 0, -4, 235, -3, 5, 0,  # $
                                  -24, 32, -41, -2, -19, -4, -3, 1,  # $
                                  0, 0, -180, 0, 0, 0, 0, 0], dtype="int")
        delay_add2 = numpy.array([0, 0, 56, 56, 40, 56, 56, 56,  # $  ;20150124
                                  56, 56, 72, 56, 56, 56, 56, 56,  # $
                                  56, 56, 56, 56, 56, 56, 56, 56,  # $
                                  568, 56, 472, -56, 56, 56, 52, 56,  # $
                                  56, 64, 56, 44, 56, 56, 56, 56,  # $
                                  56, 56, 56, 56, 56, 56, 56, 56,  # $
                                  56, 56, 56, 56, 56, 56, 56, 56,  # $
                                  56, 56, 300, 56, 56, 56, 56, 56], dtype="int")
        delay_add3 = numpy.array([-17, 0, -56, -56, -51, -56, -56, -58,  # $
                                  -57, -55, -72, -56, -58, -56, -57, -59,  # $
                                  -60, -56, -56, -57, -56, -47, -58, -55,  # $
                                  -182, -54, -159, -261, -57, 53, -56, -59,  # $
                                  -62, -65, -57, -50, -51, 656, 1824, -57,  # $
                                  -62, -57, -55, -53, -58, -58, -55, -165,  # $
                                  -57, -57, -55, -57, -62, -58, -54, -54,  # $
                                  -50, -56, -116, -56, -56, -55, -56, -56], dtype="int")

        for i in range(0, self.delay_array_hi.size):  # dr_output_antennas = 64
            self.delay_array_hi[i] += (delay_add1[i] + delay_add2[i] + delay_add3[i])

    def get_delay_value(self, sub_array=1, rdate=None):
        if sub_array == 1:
            return self.delay_array_low
        else:
            return self.delay_array_hi


