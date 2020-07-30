import os
from muserenv import *
from musersqlite import musersqlite
from muserredis import MuserRedis

class MuserObservatory(object):
    __instance = None
    def __new__(cls):
        if cls.__instance is None:
            cls.__instance = super(MuserObservatory,cls).__new__(cls)
            cls.__instance.__initialized = False
        return cls.__instance

    def __init__(self):
        if(self.__initialized): return
        self.__initialized = True
        self.__position = []

    def get_muser_observatory(self, rdate='2015-11-01 00:00:00'):
        if len(self.__position)==3:
            return self.__position[0],self.__position[1],self.__position[2]
        datasource = muserenv.get_datasource()
        if datasource == 'file':
            return self.get_from_text()
        elif datasource == 'database':
            return self.get_from_db(rdate)
        elif datasource == 'service':
            return self.get_from_service(rdate)
        # if len(self.__position)==3:
        #     return self.__position[0],self.__position[2],self.__position[2]


    def get_from_text(self):
        muser_file_path = os.path.join(muserenv.get_home_dir() + "/data", "muser.txt")
        if not os.path.exists(muser_file_path):
            return []
        file = open(muser_file_path)
        try:
            line = file.readline()
            self.__position.append(float(line.split()[0]))
            self.__position.append(float(line.split()[1]))
            self.__position.append(float(line.split()[2]))
            return self.__position
        finally:
            file.close()

    def get_from_db(self, rdate='2015-11-01 00:00:00'):
        if musersqlite is not None:
            try:
                musersqlite.connect()
                args=['observatory', rdate]
                sql = '''select refTime, theValue from t_global where keyName=? and refTime<=? order by refTime desc limit 1''' #self, sql, data
                #sql = '''select refTime, theValue from t_global where freq=? and refTime <=? order by refTime desc limit 1''' #self, sql, data

                lines = musersqlite.query(sql,args)
                if len(lines)>0:
                    line = (lines[0][1]).split(' ')
                    self.__position.append(float(line[0].encode('ascii')))
                    self.__position.append(float(line[1].encode('ascii')))
                    self.__position.append(float(line[2].encode('ascii')))
                return self.__position
            finally:
                musersqlite.close()

    def get_from_service(self, rdate='2015-11-01 00:00:00'):
        objDict = MuserRedis().muser_get_global("observatory",rdate)
        if objDict["error"] is None :
            line = objDict["result"][0]
            return float(line.split()[0]), float(line.split()[1]), float(line.split()[2])

muserobservatory = MuserObservatory()