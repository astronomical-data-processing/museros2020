import pyfits
import os, sys
class MuserFits:

    def create_fits(self, data, object, obs_date, obs_time, imagetype):
        self.hdunew = pyfits.PrimaryHDU(data=data)
        self.hdunew.header.set('OBJECT','SUN')
        self.hdunew.header.set('TELESCOP', object)
        #self.hdunew.header.set('OBS-MODE', loopMode) # observing mode: Loop / Unloop
        self.hdunew.header.set('ORIGIN', 'Mingantu Obs, NAOC of CAS')
        self.hdunew.header.set('DATA-TYP', imagetype) # "Cleaned_map")
        self.hdunew.header.set('COMMENT','FITS File of MUSER.')
        self.hdunew.header.set('OBS-DATE', obs_date, 'obs-date (CST)')
        self.hdunew.header.set('OBS-TIME', obs_time, 'obs-time (CST)')

    def append_header(self, fieldname, fieldcontent, comment):
        self.hdunew.header.set(fieldname, fieldcontent, comment)

    def append_common_header(self, freq, stocks, ra, dec, solp):

        self.hdunew.header.set('OBS-FREQ', freq, 'frequency')
        self.hdunew.header.set('STOCKS', stocks, 'polarization')

        self.hdunew.header.set('RA', ra, 'right ascension (hour)')
        self.hdunew.header.set('DEC', dec, 'declination (degree)')

        self.hdunew.header.set('SOLR', 974.65, 'optical solar radius (arcsecond)')
        self.hdunew.header.set('SOLP', solp, 'solar polar angle (degree)')
        # self.hdunew.header.set('SOLB', 0.)  # solar b0 (degree)

        self.hdunew.header.set('BUNIT',   'K ', 'disk = 10000 K')
        self.hdunew.header.set('DSKBR', 0., 'brightness of the dirty disk')

        self.hdunew.header.set('GAINE', 0.1, 'CLEAN loop gain')
        self.hdunew.header.set('TRIM', 0.1, 'CLEAN trim level')
        #self.hdunew.header.set('HIERARCH low-pass filter', 1.) # low-pass filter at 160d
        # self.hdunew.header.set('CRITER', 1000.) # actual clean criterion
        # self.hdunew.header.set('NCOMPO', 1000.) # number of clean components

        # self.hdunew.header.set('STRT-OBS', data_time) # STRT-OBS= '23:53:00.645'
        # self.hdunew.header.set('END-OBS', data_time)   # END-OBS = '23:53:10.645'
        # self.hdunew.header.set('JSTDATE', data_time)   # JSTDATE = '2014-12-01'
        # self.hdunew.header.set('JSTTIME', data_time)   # JSTTIME = '08:53:05.645'
        # self.hdunew.header.set('JST-STRT', data_time)  # JST-STRT= '08:53:00.645'
        # self.hdunew.header.set('JST-END', data_time)   # JST-END = '08:53:10.645'
        # self.hdunew.header.set('STARTFRM', 4091)
        # self.hdunew.header.set('ENDFRM', 4100)
        # self.hdunew.header.set('FRM-STAT', "1-sec obs")

        # self.hdunew.header.set('DATA-TYP', "r+1")
        # self.hdunew.header.set('ATT-10DB', "00dB")

        # self.hdunew.header.set('HOURA', 0.) # hour angle (second)
        # self.hdunew.header.set('AZIMUTH', 0.) # azimuth (degree)
        # self.hdunew.header.set('ALTITUDE', 0.) # altitude (degree)
        # self.hdunew.header.set('ZANGLE', 0.) # zenithangle (degree)

        # self.hdunew.header.set('PMAT1', 0.) # projection matrix
        # self.hdunew.header.set('PMAT2', 0.) # projection matrix
        # self.hdunew.header.set('PMAT3', 0.) # projection matrix
        # self.hdunew.header.set('PMAT4', 0.) # projection matrix

        # self.hdunew.header.set('CRVAL1',   0.00) # arcsec
        # self.hdunew.header.set('CRVAL2',   0.00) # arcsec
        # self.hdunew.header.set('CRPIX1',   256.50)
        # self.hdunew.header.set('CRPIX2',   256.50)
        # self.hdunew.header.set('CDELT1',   4.91106) # arcsec
        # self.hdunew.header.set('CDELT2',   4.91106) # arcsec
        # self.hdunew.header.set('CTYPE1',   'SOLAR-WEST')
        # self.hdunew.header.set('CTYPE2',   'SOLAR-NORTH')

        # self.hdunew.header.set('SOLR-FAC', 1.0) # radius correction factor
        # self.hdunew.header.set('NFRCAL', 1.) # number of calibration frames
        # self.hdunew.header.set('CRINPUT', 3000.00) # lower limit of clean criterion  Iterations
        # self.hdunew.header.set('CRFACTOR', 0.00) # factor clean criterion
        # self.hdunew.header.set('MBEAMC', "YES") # main beam correction
        # self.hdunew.header.set('DISKRSTR', "YES") # disk restoration

        # self.hdunew.header.set('DDOFF1', 0.) # x-offset of the dirty disk
        # self.hdunew.header.set('DDOFF2', 0.) # y-offset of the dirty disk
        # self.hdunew.header.set('DDCORR', 0.) # correlation between dirty disk and model


    def write_fits(self, outdir, filename): #, data, ra, dec, date_time):

        hdunewlistnew = pyfits.HDUList([self.hdunew])

        pathPrefix = outdir
        if pathPrefix[-1:] == '/':
            pathPrefix = pathPrefix[:-1]
        if not os.path.exists(pathPrefix):
            os.makedirs(pathPrefix)

        if os.path.exists(filename):
            os.remove(filename)
        hdunewlistnew.writeto(filename,clobber=True)
