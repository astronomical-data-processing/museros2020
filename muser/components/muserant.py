from muserenv import *
import threading

class MuserAntenna(object):
    __instance = None
    def __new__(cls):
        if cls.__instance is None:
            cls.__instance = super(MuserAntenna,cls).__new__(cls)
            cls.__instance.__initialized = False
        return cls.__instance

    def __init__(self, sub_array = 1):
        if(self.__initialized): return
        self.antennas=[40,60]
        self.muser_ant_1 = []
        self.muser_ant_2 = []
        self.antenna1_loaded = False
        self.antenna2_loaded = False
        self.lock = threading.Lock()
        self.sub_array = sub_array
        #self.env = MuserEnv()
        self.Flag_Ant = { '2015-11-01':[8, 9,10, 11, 12, 13, 19, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 34, 35, 36, 37, 38, 39],
                 '2014-11-11':[4, 7, 10, 11,12,13,16,17,18, 19, 24, 25, 26, 36, 38, 39]}

    def get_muser_antennas_position(self, sub_array, rdate):
        if muserenv.get_datasource() == 'file':
            if sub_array == 1:
                if self.antenna1_loaded == True:
                    return self.muser_ant_1
            elif sub_array == 2:
                if self.antenna2_loaded == True:
                    return self.muser_ant_2

            # get from file
            antenna_file_path = os.path.join(muserenv.get_home_dir() + "/data", "ant_pos_%1d.txt" % (sub_array))
            file = open(antenna_file_path)
            self.lock.acquire()
            try:
                self.muser_ant_1 = []
                for i in range(self.antennas[sub_array -1]):
                    line = file.readline()
                    if sub_array ==1:
                        self.muser_ant_1.append(float(line.split()[0]))
                        self.muser_ant_1.append(float(line.split()[1]))
                        self.muser_ant_1.append(float(line.split()[2]))
                    elif sub_array ==2:
                        self.muser_ant_2.append(float(line.split()[0]))
                        self.muser_ant_2.append(float(line.split()[1]))
                        self.muser_ant_2.append(float(line.split()[2]))
                if sub_array ==1:
                    self.antenna1_loaded = True
                    return self.muser_ant_1
                elif sub_array ==2:
                    self.antenna2_loaded = True
                    return self.muser_ant_2
            finally:
                file.close()
                self.lock.release()

        elif  muserenv.get_datasource() == 'database':
            try:
                musersqlite.connect()
                sql = '''select refTime, theValue from p_antenna_position where freq=? and refTime<=? order by refTime desc limit 1'''
                arguments=[sub_array,rdate]
                antenna = musersqlite.query(sql, arguments)
                # print "Antenna:", antenna
                if len(antenna)>0:
                    str_ant = antenna[0][1].split(',')
                    muser_ant = []
                    for a in range(len(str_ant)):
                        muser_ant.append(float(str_ant[a].encode('ascii')))
                    return muser_ant
                else:
                    return []
            finally:
                musersqlite.close()

        elif muserenv.get_datasource() == 'service':
            objDict = MuserRedis().muser_get_position(sub_array,rdate)
            if objDict["error"] is None :
                str_ant = objDict["result"][0].split(',')
                muser_ant = []
                for a in range(len(str_ant)):
                    muser_ant.append(float(str_ant[a].encode('ascii')))
                return muser_ant
            else:
                return []


    def get_flag_antenna(self, sub_array, rdate):
        if self.Flag_Ant.has_key(rdate):
            return self.Flag_Ant[rdate]

        if muserenv.get_datasource() == 'database':
            if musersqlite is not None:
                try:
                    musersqlite.connect()
                    args=([sub_array,rdate])
                    sql = '''select refTime, theValue from p_antenna_flag where freq=? and refTime <=? order by refTime desc limit 1''' #self, sql, data

                    lines = musersqlite.query(sql,args)
                    if len(lines)>0:
                        str_ant = lines[0][1].split(',')
                        #print "LINE:",line
                        muser_ant = []
                        for a in range(len(str_ant)):
                            muser_ant.append(float(str_ant[a].encode('ascii')))
                        return muser_ant
                finally:
                    musersqlite.close()
        elif muserenv.get_datasource() == 'service':
            objDict = MuserRedis().muser_get_flag(sub_array,rdate)
            if objDict["error"] is None :
                str_ant = objDict["result"][0].split(',')
                muser_ant = []
                for a in range(len(str_ant)):
                    muser_ant.append(float(str_ant[a].encode('ascii')))
                return muser_ant
            else:
                return []

        # elif source==1:
        #     musersqlite


muserantenna = MuserAntenna()