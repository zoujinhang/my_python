import astropy.time as time
#import astropy.units as u
#import numpy as np

class GBMtime(object):
	def __init__(self):
		self.utc_start_time = '2008-08-07T03:35:44.0'
		self.mjd_start_time = 51910

	@classmethod   #加上这个就可以使下面的过程直接运行。

	def met_to_utc(self,met):
		if (met <= 252460801.000):
			utc_tt_diff = 65.184
		elif (met <= 362793602.000):
			utc_tt_diff = 66.184
		elif (met <= 457401603.000):
			utc_tt_diff = 67.184
		elif (met <=  504921604.000):
			utc_tt_diff = 68.184
		else:
			utc_tt_diff = 69.184

		mjdutc = ((met - utc_tt_diff) / 86400.0) + 51910 + 0.0007428703703
		met1 = time.Time(mjdutc,scale= 'utc',format = 'mjd')

		return met1

	@classmethod

	def utc_to_met(self,utc0):
		tt_time = time.Time(utc0, format='fits', scale='utc').mjd
		mmt = (tt_time - 0.0007428703703 - 51910) * 86400.0
		if mmt <= (252460801.000 - 65.184):
			dt = 65.184
		elif mmt <= (362793602.000 - 66.184):
			dt = 66.184
		elif mmt <= (457401603.000 - 67.184):
			dt = 67.184
		elif mmt <= (504921604.000 - 68.184):
			dt = 68.184
		else:
			dt = 69.184
		met = mmt + dt
		return met

	@classmethod

	def utc_time(self,utc0):
		tt = time.Time(utc0,format = 'fits',scale = 'utc')
		return tt

#@staticmethod  #这个是静态变量的意思。
#t = '2018-06-30T00:00:00'

#ee = GBMtime.utc_time(t)
#print(ee)








