import numpy as np
import re
import os
from astropy.time import Time
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from astropy.io import fits
from zjh_gbm_time import GBMtime
from zjh_download import download_all_in_one_path
import zzh_py3_file as zhfl
import astropy.units as u
from zjh_gbm import *


class Fermi_daily_geometry(object):
	def __init__(self,time):

		self.topdir = '/home/laojin/gbm_daily_database/'
		self.resultdir = '/home/laojin/shiyan/'
		self.time = Time(time).fits
		time_list = re.split('[T,-]',self.time)
		self.localdir = self.topdir + time_list[0] + '/' + time_list[1] + '/' + time_list[2] + '/'
		filef = 'glg_poshist_all_'+time_list[0][2:]+time_list[1]+time_list[2]
		if(os.path.exists(self.localdir) == False ):
			targetdir = '/fermi/data/gbm/daily/' + time_list[0]+'/'+time_list[1] + '/' + time_list[2] + '/current/'
			download_all_in_one_path(targetdir,self.localdir)

		filelist = zhfl.findfile(self.localdir,filef+'*')
		if((filelist == 0) | (filelist == False)):
			targetdir = '/fermi/data/gbm/daily/' + time_list[0]+'/'+time_list[1] + '/' + time_list[2] + '/current/'
			download_all_in_one_path(targetdir,self.localdir)
			filelist = zhfl.findfile(self.localdir,filef+ '*')
		self.filename = filelist[0]
		self.filelink = os.path.join(self.localdir,self.filename)
		self.met = GBMtime.utc_to_met(self.time)
	def open_fit(self):
		f = fits.open(self.filelink)
		time = f[1].data.field(0)
		qsj1 = f[1].data.field(1)
		qsj2 = f[1].data.field(2)
		qsj3 = f[1].data.field(3)
		qsj4 = f[1].data.field(4)
		pos_x = f[1].data.field(8)
		pos_y = f[1].data.field(9)
		pos_z = f[1].data.field(10)
		sc_lat = f[1].data.field(14)
		sc_lon = f[1].data.field(15)
		return time,qsj1,qsj2,qsj3,qsj4,pos_x,pos_y,pos_z,sc_lat,sc_lon
	def find_right_list(self):
		time,qsj1,qsj2,qsj3,qsj4,pos_x,pos_y,pos_z,sc_lat,sc_lon = self.open_fit()
		t = (time - self.met)**2
		t = np.array(t)
		index = np.where(t == np.min(t))
		print(t[index][0])
		print('met: ',time[index][0])
		qsj = np.array([qsj1[index][0],qsj2[index][0],qsj3[index][0],qsj4[index][0]])
		pos = np.array([pos_x[index][0],pos_y[index][0],pos_z[index][0]])
		sc = np.array([sc_lat[index][0],sc_lon[index][0]])
		return qsj,pos,sc

	def get_gbmgeometry(self,plot = True):
		qsj,pos,sc = self.find_right_list()
		myGBM = GBM(qsj,sc_pos=pos*u.m,gbm_time=self.time)
		if plot:
			myGBM.detector_plot(radius = 10,lat_0 = 0,lon_0 = 180,show_bodies = True,BGO= False)
			resultname = 'gbm_'+self.time + '.png'
			plt.title(time)
			print(self.resultdir+resultname)
			plt.savefig(self.resultdir+resultname)
		return myGBM

	def plot_position_on_earth(self,plot_fermi = True):
		time1,qsj1,qsj2,qsj3,qsj4,pos_x,pos_y,pos_z,sc_lat,sc_lon = self.open_fit()
		t = (time1 - self.met)**2
		t = np.array(t)
		index = np.where(t == np.min(t))
		sc = np.array([sc_lat[index][0],sc_lon[index][0]])
		lat_0 = np.array(
			[-30.0, -22.6, 2.5, 5.2, 5.2, 4.6, 0.7, -8.6, -9.9, -12.5, -21.7, -30.0, -52.5, -52.5, -45.9, -37.5,
			 -35.0, -30])
		lon_0 = np.array(
			[33.9, 24.5, -18.6, -25.7, -36.0, -42.0, -58.8, -93.1, -97.5, -98.5, -92.1, -86.1, -45.0, -30.0, 0.0,
			 30.0, 37.5, 33.9])
		lon_1 = np.array([33.900000,24.306791,20.475922,16.612362,12.409036,2.1819403,-9.5971985,-22.151855,-36.133575,-40.414093,-60.145325,-65.28494,-76.805084,-89.97073,-94.29413,-94.20526,-91.904236,-90.0206,-89.500244,-86.100000, -45.0,-30.0, 0.0,30.0, 37.5, 33.9])
		lat_1= np.array([-30.000000,-25.623365,-23.76377,-22.037643,-19.97555,-15.081318,-9.626598,-3.8377674,1.486362,1.8389976,-0.52080667,-1.0369452,-4.446427,-9.900604,-14.818653,-18.087845,-21.907251,-24.687246,-25.313162,-30.000000, -52.5, -52.5, -45.9, -37.5,-35.0, -30])

		fig = plt.figure(figsize=(20, 10))
		ax = fig.add_subplot(111)
		mp = Basemap(projection='moll', lon_0=0, lat_0=0, resolution='l',ax = ax)
		mp.drawcoastlines()
		mp.drawcountries()
		mp.drawmapboundary()
		mp.drawmeridians(np.arange(0, 360, 30))
		mp.drawparallels(np.arange(-90, 90, 30))
		x, y = mp(lon_0, lat_0)
		mp.plot(x, y,linewidth = 5 ,color='#aa2116')
		x,y = mp(lon_1,lat_1)
		mp.plot(x,y,linewidth = 5,color ='#5c7a29')
		point_x,point_y = mp(-50,-30)
		plt.text(point_x,point_y,'SAA',color = '#892f1b',size = 50)
		lon, lat = mp(sc_lon, sc_lat)
		#mp.plot(lon, lat, ',', color='#f47920')
		mp.plot(lon, lat, '.', color='#905d1d',markersize = 2)
		if plot_fermi:
			lon1, lat1 = mp(sc[1], sc[0])
			mp.plot(lon1, lat1, '.', color='#f58220', markersize=30)
			plt.text(lon1, lat1, '  Fermi', color='#f58220', size=30)
		return mp




















