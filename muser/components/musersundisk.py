import sys, os, datetime, time, math
from muserframe import *
from muserenv import *
import numpy as np

logger = logging.getLogger('muser')


class MuserSunDisk(object):
    __instance = None
    def __new__(cls):
        if cls.__instance is None:
            cls.__instance = super(MuserSunDisk,cls).__new__(cls)
            cls.__instance.__initialized = False
        return cls.__instance

    def __init__(self, cached=False):
        if (self.__initialized): return
        self.__initialized = True
        self.cached = cached
        if cached == True:
            self.cachedata()


    def cachedata(self):
        sundisk_file_path = os.path.join(muserenv.get_home_dir() + "/data", "quietsundisk.txt")
        # sundisk_file_path ="/home/fengyong/data/quietsundisk2400_4375MHz.txt"
        sundisk_file = open(sundisk_file_path, "r")
        self.sundisk_dic = {}
        try:
            while True:
                linekey = sundisk_file.readline()
                linevalue = sundisk_file.readline()
                if linekey and linevalue:
                    linekey = float(linekey.strip())
                    if linekey == frequency:
                        value = map(eval, linevalue.strip().split(',')[1:])
                        npvalue = []
                        for i in range(200):
                            npvalue.append(np.float32(value[i]))
                        # npvalue = npvalue/npvalue[0]
                        self.sundisk_dic[linekey] = np.array(npvalue)
                else:
                    break
                    # print type(self.sundisk_dic[1725.00][0]), max(self.sundisk_dic[1725.00]), min(self.sundisk_dic[1725.00])
        finally:
            sundisk_file.close()

    def getdata(self, frequency):
        if self.cached== True:
            return self.sundisk_dic[frequency]
        else:
            sundisk_file_path = os.path.join(muserenv.get_home_dir() + "/data", "quietsundisk.txt")
            # sundisk_file_path ="/home/fengyong/data/quietsundisk2400_4375MHz.txt"
            sundisk_file = open(sundisk_file_path, "r")
            self.sundisk_dic = {}
            try:
                while True:
                    linekey = sundisk_file.readline()
                    linevalue = sundisk_file.readline()
                    if linekey and linevalue:
                        linekey = float(linekey.strip())
                        if linekey == frequency:
                            value = map(eval, linevalue.strip().split(',')[1:])
                            npvalue = []
                            for i in range(200):
                                if i<100:
                                    npvalue.append(np.float32(1.))
                                else:
                                    npvalue.append(np.float32(0.))
                                # npvalue.append(np.float32(value[i]))
                            nvalue = np.array(npvalue)
                            max = np.max(nvalue)
                            nvalue = nvalue/max
                            return nvalue
                    else:
                        break
                        # print type(self.sundisk_dic[1725.00][0]), max(self.sundisk_dic[1725.00]), min(self.sundisk_dic[1725.00])
            finally:
                sundisk_file.close()

if __name__ == '__main__':
    sundisk = MuserSunDisk()
    sundisk.readdata()
    frequency = 400.00
    print sundisk.sundisk_dic
    # num = len(sundisk.sundisk_dic[frequency])
    #
    # sundisk.draw_sundiskmodel(frequency,num)
