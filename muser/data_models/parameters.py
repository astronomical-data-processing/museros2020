"""We use the standard kwargs mechanism for arguments. For example::

    kernelname = get_parameter(kwargs, "kernel", "2d")
    oversampling = get_parameter(kwargs, "oversampling", 8)
    padding = get_parameter(kwargs, "padding", 2)

The kwargs may need to be passed down to called functions.

All functions possess an API which is always of the form::

      def processing_function(idatastruct1, idatastruct2, ..., *kwargs):
         return odatastruct1, odatastruct2,... other

Processing parameters are passed via the standard Python kwargs approach.

Inside a function, the values are retrieved can be accessed directly from the
kwargs dictionary, or if a default is needed a function can be used::

    log = get_parameter(kwargs, 'log', None)
    vis = get_parameter(kwargs, 'visibility', None)

Function parameters should obey a consistent naming convention:

=======  =======
Name     Meaning
=======  =======
vis      Name of Visibility
sc       Name of Skycomponent
gt       Name of GainTable
conf     Name of Configuration
im       Name of input image
qa       Name of quality assessment
log      Name of processing log
=======  =======

If a function argument has a better, more descriptive name e.g. normalised_gt, newphasecentre, use it.

Keyword=value pairs should have descriptive names. The names should be lower case with underscores to separate words:

====================    ==================================  ========================================================
Name                    Meaning                             Example
====================    ==================================  ========================================================
loop_gain               Clean loop gain                     0.1
niter                   Number of iterations                10000
eps                     Fractional tolerance                1e-6
threshold               Absolute threshold                  0.001
fractional_threshold    Threshold as fraction of e.g. peak  0.1
G_solution_interval     Solution interval for G term        100
phaseonly               Do phase-only solutions             True
phasecentre             Phase centre (usually as SkyCoord)  SkyCoord("-1.0d", "37.0d", frame='icrs', equinox='J2000')
spectral_mode           Visibility processing mode          'mfs' or 'channel'
====================    ==================================  ========================================================

"""

__all__ = ['muser_path', 'muser_data_path', 'get_parameter', 'IF_BANDWIDTH', 'LOOP_MODE_LOW', 'NON_LOOP_MODE_LOW',
           'LOOP_MODE_HIGH','NON_LOOP_MODE_HIGH']

import logging
import os

log = logging.getLogger('logger')

# Dictionary of Bandwidth of Muser
IF_BANDWIDTH = {0x33333333: 25000000, 0x77777777: 12500000, 0x88888888: 6250000, 0xbbbbbbbb: 3125000,
                0xcccccccc: 1562500}

# Dictionary of Bandwidth of Muser-I
LOOP_MODE_LOW = {0x0: (0, 0, 400000000, 1), 0x1: (1, 0, 400000000, 1), 0x2: (0, 16, 800000000, 2),
                 0x3: (1, 16, 800000000, 2), 0x4: (0, 32, 1200000000, 3), 0x5: (1, 32, 1200000000, 3),
                 0x6: (0, 48, 1600000000, 4), 0x7: (1, 48, 1600000000, 4)}

# Dictionary of Non-loop mode of Muser-I
NON_LOOP_MODE_LOW = {0x3333: (0, 400000000, 1), 0x7777: (16, 800000000, 2),
                     0xbbbb: (32, 1200000000, 3), 0xcccc: (48, 1600000000, 4)}

# Dictionary of Loop mode of Muser-II
LOOP_MODE_HIGH = {0X0: (0, 0, 2000000000), 0X1: (1, 0, 2000000000), 0X2: (0, 16, 2400000000),
                  0X3: (1, 16, 2400000000), 0X4: (0, 32, 2800000000), 0X5: (1, 32, 2800000000),
                  0X6: (0, 48, 3200000000), 0X7: (1, 48, 3200000000), 0X8: (0, 64, 3600000000),
                  0X9: (1, 64, 3600000000), 0XA: (0, 80, 4000000000), 0XB: (1, 80, 4000000000),
                  0XC: (0, 96, 4400000000), 0XD: (1, 96, 4400000000), 0XE: (0, 112, 4800000000),
                  0XF: (1, 112, 4800000000), 0X10: (0, 128, 5200000000), 0X11: (1, 128, 5200000000),
                  0X12: (0, 144, 5600000000), 0X13: (1, 144, 5600000000), 0X14: (0, 160, 6000000000),
                  0X15: (1, 160, 6000000000), 0X16: (0, 176, 6400000000), 0X17: (1, 176, 6400000000),
                  0X18: (0, 192, 6800000000), 0X19: (1, 192, 6800000000), 0X1A: (0, 208, 7200000000),
                  0X1B: (1, 208, 7200000000), 0X1C: (0, 224, 7600000000), 0X1D: (1, 224, 7600000000),
                  0X1E: (0, 240, 8000000000), 0X1F: (1, 240, 8000000000), 0X20: (0, 256, 8400000000),
                  0X21: (1, 256, 8400000000), 0X22: (0, 272, 8800000000), 0X23: (1, 272, 8800000000),
                  0X24: (0, 288, 9200000000), 0X25: (1, 288, 9200000000), 0X26: (0, 304, 9600000000),
                  0X27: (1, 304, 9600000000), 0X28: (0, 320, 10000000000), 0X29: (1, 320, 10000000000),
                  0X2A: (0, 336, 10400000000), 0X2B: (1, 336, 10400000000), 0X2C: (0, 352, 10800000000),
                  0X2D: (1, 352, 10800000000), 0X2E: (0, 368, 11200000000), 0X2F: (1, 368, 11200000000),
                  0X30: (0, 384, 11600000000), 0X31: (1, 384, 11600000000), 0X32: (0, 400, 12000000000),
                  0X33: (1, 400, 12000000000), 0X34: (0, 416, 12400000000), 0X35: (1, 416, 12400000000),
                  0X36: (0, 432, 12800000000), 0X37: (1, 432, 12800000000), 0X38: (0, 448, 13200000000),
                  0X39: (1, 448, 13200000000), 0X3A: (0, 464, 13600000000), 0X3B: (1, 464, 13600000000),
                  0X3C: (0, 480, 14000000000), 0X3D: (1, 480, 14000000000), 0X3E: (0, 496, 14400000000),
                  0X3F: (1, 496, 14400000000), 0X40: (0, 512, 14600000000), 0X41: (1, 512, 14600000000)}
NON_LOOP_MODE_HIGH = {0x0: (0, 2000000000, 1), 0x101: (16, 2400000000, 2), 0x202: (32, 2800000000, 3),
                      0x303: (48, 3200000000, 4), 0x404: (64, 3600000000, 5), 0x505: (80, 4000000000, 6),
                      0x606: (96, 4400000000, 7), 0x707: (112, 4800000000, 8), 0x808: (128, 5200000000, 9),
                      0x909: (144, 5600000000, 10), 0xa0a: (160, 6000000000, 11), 0xb0b: (176, 6400000000, 12),
                      0xc0c: (192, 6800000000, 13), 0xd0d: (208, 7200000000, 14), 0xe0e: (224, 7600000000, 15),
                      0xf0f: (240, 8000000000, 16), 0x1010: (256, 8400000000, 17), 0x1111: (272, 8800000000, 18),
                      0x1212: (288, 9200000000, 19), 0x1313: (304, 9600000000, 20), 0x1414: (320, 10000000000, 21),
                      0x1515: (336, 10400000000, 22), 0x1616: (352, 10800000000, 23),
                      0x1717: (368, 11200000000, 24),
                      0x1818: (384, 11600000000, 25), 0x1919: (400, 12000000000, 26),
                      0x1a1a: (416, 12400000000, 27),
                      0x1b1b: (432, 12800000000, 28), 0x1c1c: (448, 13200000000, 29),
                      0x1d1d: (464, 13600000000, 30),
                      0x1e1e: (480, 14000000000, 31), 0x1f1f: (496, 14400000000, 32),
                      0x2020: (512, 14600000000, 33)}


def muser_path(path):
    """Converts a path that might be relative to muser root into an
    absolute path::

        muser_data_path('models/SKA1_LOW_beam.fits')
        '/Users/wangfeng/work/muser2020/data/models/SKA1_LOW_beam.fits'

    :param path:
    :return: absolute path
    """
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)) + "/../../")
    muserhome = os.getenv('MUSER', project_root)
    return os.path.join(muserhome, path)


def muser_data_path(path=None, check=True):
    """Converts a path that might be relative to the muser data directory into an
    absolute path::

        muser_data_path('models/SKA1_LOW_beam.fits')
        '/Users/wangfeng/work/muser2020/data/models/SKA1_LOW_beam.fits'

    The data path default is muser_path('data') but may be overriden with the environment variable muser_DATA.

    :param path:
    :return: absolute path
    """
    muser_data_home = os.getenv('MUSER_DATA', None)
    if muser_data_home is None:
        dp = muser_path('configurations/')
    else:
        dp = muser_data_home
    if check:
        if not os.path.exists(dp):
            raise EnvironmentError("MUSER data directory {} does not exist".format(dp))
    dp = os.path.join(dp, path)
    return dp

def get_parameter(kwargs, key, default=None):
    """ Get a specified named value for this (calling) function

    The parameter is searched for in kwargs

    :param kwargs: Parameter dictionary
    :param key: Key e.g. 'loop_gain'
    :param default: Default value
    :return: result
    """

    if kwargs is None:
        return default

    value = default
    if key in kwargs.keys():
        value = kwargs[key]
    return value
