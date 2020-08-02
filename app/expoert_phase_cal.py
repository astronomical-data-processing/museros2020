"""Export 

"""

import os, sys

path = os.path.abspath(os.path.dirname(__file__))

from argparse import *
import datetime
import traceback
import logging
from muser.data_models.muser_data import MuserData
from muser.data_models.parameters import muser_path, muser_data_path
from rascil.processing_components.visibility.coalesce import convert_visibility_to_blockvisibility, \
    convert_blockvisibility_to_visibility


log = logging.getLogger('muser')

class Phase:
    def __init__(self, sub_array=None, is_loop_mode=None, obs_date_time=None, frame_number=None, calibration_sequence=None, filename=None, debug=0):
        self.sub_array = sub_array
        self.frame_number = frame_number
        self.calibration_sequence = calibration_sequence
        self.data_source = 0
        self.obs_date_time = None
        self.debug = debug
        self.is_loop_mode = is_loop_mode
        if obs_date_time is not None and len(obs_date_time)>0:
            log.info("Searched phase data date and time: %s" % obs_date_time)
            self.obs_date_time = valid_date(obs_date_time)
            self.year = self.obs_date_time.date().year
            self.month = self.obs_date_time.date().month
            self.day = self.obs_date_time.date().day
            self.hour = self.obs_date_time.time().hour
            self.minute = self.obs_date_time.time().minute

        self.filename = None
        self.filename = filename.strip()

    def Calibration(self):
        if self.frame_number < 0:
            log.error("You should input a positive number!")
            return False

        muser_calibration = MuserData(self.sub_array)
        if self.obs_date_time is not None:
            muser_calibration.set_data_date_time(self.year, self.month, self.day, self.hour, self.minute,0,0,0,0)
        if len(self.filename.strip())>0:
            muser_calibration.input_file_name=self.filename
        else:
            muser_calibration.input_file_name=''
        if self.debug == 1:
            log.info('Reading Visibility Data of calibration......')
        if muser_calibration.open_data_file() == False:
            print 'Error: cannot find the data file.'
            exit


        muser_calibration.skip_frames(self.frame_number)
        self.last_sub_band = -1
        self.last_polarization = -1

        if self.is_loop_mode == True:
            frame_NUM = muser_calibration.frame_number * 2

            calibration_Data = np.ndarray(
                shape=(muser_calibration.frame_number, muser_calibration.polarization_number,
                       muser_calibration.antennas * (muser_calibration.antennas - 1) // 2, 16),
                dtype=complex)
        else:
            frame_NUM = 1
            print
            calibration_Data = np.ndarray(
                shape=(muser_calibration.antennas * (muser_calibration.antennas - 1) // 2, 16),
                dtype=complex)

        for i in range(frame_NUM):
            if (muser_calibration.read_one_frame() == False):  # 32*8bits
                log.error('Cannot read a frame.')
                return False
            muser_calibration.read_data()

            self.year = muser_calibration.current_frame_time.year
            self.month = muser_calibration.current_frame_time.month
            self.day = muser_calibration.current_frame_time.day

            if self.debug == 1:
                log.info("Reading No. %d %s %d %d" % (i, muser_calibration.current_frame_time.get_string(), muser_calibration.sub_band, muser_calibration.polarization))
            #According to the DISCUSSION with LJ.Chen
            #muser_calibration.delay_process('satellite')
            # Delay processing for satellite
            if self.sub_array == 2:
                if muser_calibration.current_frame_header.strip_switch == 0xCCCCCCCC:
                    muser_calibration.delay_process("satellite")

            self.last_sub_band =  muser_calibration.sub_band
            self.last_polarization = muser_calibration.polarization

            bl = 0
            for antenna1 in range(0, muser_calibration.antennas - 1):
                for antenna2 in range(antenna1 + 1, muser_calibration.antennas):
                    for channel in range(0, muser_calibration.sub_channels):
                        if self.is_loop_mode == True:
                            calibration_Data[muser_calibration.sub_band][muser_calibration.polarization][bl][channel] = muser_calibration.baseline_data[bl][channel]
                        else:
                            calibration_Data[bl][channel] = muser_calibration.baseline_data[bl][channel]
                    bl = bl + 1

        file_name = self.env.cal_file(self.sub_array,self.year, self.month, self.day, self.calibration_sequence)
        if self.debug == 1:
            log.info("Writing to file: " + os.path.basename(file_name))
        calibration_Data.tofile(file_name)
        if self.debug == 1:
            log.info("Exportphase done.")
        return True

def valid_date(s):
    try:
        s = s.strip()
        split_s = string.split(s, ' ')
        if len(split_s) == 1:
            return datetime.datetime.strptime(s, "%Y-%m-%d")
        elif len(split_s) == 2:
            return datetime.datetime.strptime(s, "%Y-%m-%d %H:%M:%S")
        elif len(split_s) == 3:
            s = string.join(split_s, ' ')
            return datetime.datetime.strptime(s, "%Y-%m-%d %H:%M:%S %f")
        else:
            msg = "Not a valid date: '{0}'.".format(s)
            raise ArgumentTypeError(msg)
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise ArgumentTypeError(msg)

def exportphase (
    muser=None,
    mode=None,
    inputfile=None,
    start=None,
    frame=1,
    calibration=0,
    debug=None,
    ):

    try:
        log.origin('exportphase')
        # -----------------------------------------'
        # beginning of importmiriad implementation
        # -----------------------------------------
        #obsFileName = valid_date(start)

        cal = Phase(muser, mode, start, frame, calibration, inputfile, debug)
        cal.Calibration()
    except Exception, e:
        print traceback.format_exc()
        raise
    	log.post("Failed to export muser phase file", "ERROR")
    return


