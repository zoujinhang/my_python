import re
import os
import zzh_py3_file as zhfl
from zjh_download import download_all_in_one_path
#from ftplib import FTP
from astropy.io import fits
from astropy.time import Time
#import numpy as np
#from gbmgeometry import *
from zjh_gbm import *
import astropy.units as u
import matplotlib.pyplot as plt
#from astropy.coordinates import SkyCoord


def get_gbmgeometry(time,plot = True):
	topdir = '/home/laojin/gbm_daily_database/'
	resultdir = '/home/laojin/shiyan/'
	if isinstance(time,str) == False:
		print('E:the type you input is not str.')
		return False
	elif(len(time.split('T')) != 2):
		print('E:the format of time should like this:'+'\'2001-01-01T00:00:00\'')
		return False

	else:
		time_list = re.split('[T,-]',time)
		localdir = topdir + time_list[0] + '/' + time_list[1] + '/' + time_list[2] + '/'
		filef = 'glg_poshist_all_'+time_list[0][2:]+time_list[1]+time_list[2]
		if(os.path.exists(localdir) == False ):
			targetdir = '/fermi/data/gbm/daily/' + time_list[0]+'/'+time_list[1] + '/' + time_list[2] + '/current/'
			download_all_in_one_path(targetdir,localdir)
			#os.makedirs(localdir)
			#dowelocad_fun(time_list,localdir)

		filelist = zhfl.findfile(localdir,filef+'*')
		if((filelist == 0) | (filelist == False)):
			targetdir = '/fermi/data/gbm/daily/' + time_list[0]+'/'+time_list[1] + '/' + time_list[2] + '/current/'
			download_all_in_one_path(targetdir,localdir)
			#dowelocad_fun(time_list,localdir)
			#filelist = zhfl.findfile(localdir,filef+ '*')
		filename = filelist[0]
		filelink = os.path.join(localdir,filename)
		met = get_met_from_utc(time)
		qsj,pos,sc = find_right_list(filelink,met)
		myGBM = GBM(qsj,sc_pos=pos*u.m,gbm_time=time)
		if plot:
			myGBM.detector_plot(radius = 10,lat_0 = 0,lon_0 = 180,show_bodies = True,BGO= False)
			resultname = 'gbm_'+time + '.png'
			plt.title(time)
			print(resultdir+resultname)
			plt.savefig(resultdir+resultname)
		return myGBM

def dowelocad_fun(time_list,localdir):
	star = os.getcwd()
	target = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'+time_list[0]+'/'+time_list[1] + '/' + time_list[2] + '/current/' + 'glg_poshist_all_'+time_list[0][2:]+time_list[1]+time_list[2]+'_v00.fit'
	#print(target)
	os.chdir(localdir)
	t_link = 'wget --quiet --show-progress --read-timeout=5 --tries=0 ' + target
	os.system(t_link)
	os.chdir(star)

def get_met_from_utc(time):
	tt_time = Time(time,format = 'fits',scale = 'utc').mjd
	mmt = (tt_time-0.0007428703703-51910)*86400.0
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

def open_fit(file_link):
	f = fits.open(file_link)
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

def find_right_list(file_link,met):
	time,qsj1,qsj2,qsj3,qsj4,pos_x,pos_y,pos_z,sc_lat,sc_lon = open_fit(file_link)
	t = (time - met)**2
	t = np.array(t)
	index = np.where(t == np.min(t))
	print(t[index][0])
	print('met: ',time[index][0])
	qsj = np.array([qsj1[index][0],qsj2[index][0],qsj3[index][0],qsj4[index][0]])
	pos = np.array([pos_x[index][0],pos_y[index][0],pos_z[index][0]])
	sc = np.array([sc_lat[index][0],sc_lon[index][0]])
	return qsj,pos,sc

#t = '2018-01-03T08:20:25'
#t = '2017-03-24T07:53:55.36'
#myGBM = get_gbmgeometry(t,plot = False)
#grb = SkyCoord(ra = 217.43 , dec= -62.68,frame='icrs',unit='deg')
#separation = myGBM.get_separation(grb)
#separationlist = np.array([i for i in separation['Separation']])
#nn = separationlist <= 60
#print('----------------------------')
#print(nn.astype(int))



