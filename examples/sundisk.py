import numpy

def read_sun_disk_data():
    sundisk_file_path = "quietsundisk400_2000MHz.txt"
    sundisk_file = open(sundisk_file_path, "r")
    sundisk_dic = {}
    try:
        while True:
            linekey = sundisk_file.readline()
            linevalue = sundisk_file.readline()
            if linekey and linevalue:
                linekey = float(linekey.strip())
                value = list(map(eval, linevalue.strip().split(',')[1:]))
                npvalue = []
                for i in range(200):
                    npvalue.append(numpy.float32(value[i]))
                sundisk_dic[linekey] = numpy.array(npvalue)
            else:
                return sundisk_dic
    finally:
        sundisk_file.close()


if __name__ == '__main__':
    sun_disk = read_sun_disk_data()
    print(sun_disk)
    print(sun_disk[1975.0])