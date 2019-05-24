import astropy.time as time
import re
import numpy as np
import os
from zjh_download import download_all_in_one_path
import zzh_py3_file as zhfl
from astropy.io import fits
from zjh_gbm_time import GBMtime

def get_daily_data(time0,detector):
	time0[0] = time.Time(time0[0]).fits
	time0[1] = time.Time(time0[1]).fits
	time_list = get_time_list(time0[0],time0[1])
	path_list = get_path_list(time_list)
	filef_list = get_file_list(time_list,detector)
	metb_list = get_metb_list(time_list)
	t = np.array([])
	ch = np.array([])
	ch_n = []
	e1 = []
	e2 = []
	for index,path in enumerate(path_list):
		time1,ch1,ch_n,e1,e2 = get_data(path,filef_list[index])
		if(len(metb_list) == 1):
			t_index = np.where((time1 >= GBMtime.utc_to_met(time0[0])) & (time1 <= GBMtime.utc_to_met(time0[1])))
			time1 = time1[t_index]
			ch1 = ch1[t_index]
			t = np.concatenate([t,time1])
			ch = np.concatenate([ch,ch1])
		else:
			if(index == 0):
				t_index = np.where((time1 >= GBMtime.utc_to_met(time0[0])) & (time1 < metb_list[index + 1]))

				time1 = time1[t_index]
				ch1 = ch1[t_index]
				t = np.concatenate([t,time1])
				ch = np.concatenate([ch,ch1])
			elif(index == len(time_list)-1):
				t_index = np.where((time1 >= metb_list[index]) & (time1 <= GBMtime.utc_to_met(time0[1])))
				time1 = time1[t_index]
				ch1 = ch1[t_index]
				t = np.concatenate([t,time1])
				ch = np.concatenate([ch,ch1])
			else:
				t_index = np.where((time1 >= metb_list[index]) & (time1 < metb_list[index + 1]))
				time1 = time1[t_index]
				ch1 = ch1[t_index]
				t = np.concatenate([t,time1])
				ch = np.concatenate([ch,ch1])
	return t,ch,ch_n,e1,e2

def get_metb_list(time_list):
	met = []
	for t in time_list:
		t_list = re.split(r'[T : -]',t)
		time_c = t_list[0]+'-'+t_list[1]+'-'+t_list[2]+'T'+t_list[3]+':00:00'
		tt = GBMtime.utc_to_met(time_c)
		met.append(tt)
	return np.array(met)

def get_data(path,filef):
	topdir = '/home/laojin/gbm_daily_database/'
	localdir = topdir + path
	if(os.path.exists(localdir) == False):
		targetdir = '/fermi/data/gbm/daily/' + path + '/current/'
		download_all_in_one_path(targetdir,localdir)
	filelist = zhfl.findfile(localdir,filef+'*')
	if((filelist == 0) | (filelist == False)):
		targetdir = '/fermi/data/gbm/daily/' + path + '/current/'
		download_all_in_one_path(targetdir,localdir)
		filelist = zhfl.findfile(localdir, filef + '*')
	filename = filelist[0]
	filelink = os.path.join(localdir,filename)
	return open_daily_data_file(filelink)

def get_path_list(time_list):
	path_list = []
	for t in time_list:
		t_list = re.split(r'[T : -]',t)
		localdir = t_list[0]+'/'+t_list[1]+'/'+t_list[2]+'/'
		path_list.append(localdir)
	return path_list
def get_file_list(time_list,detector):
	file_list = []
	for f in time_list:
		f_list = re.split(r'[T : -]',f)
		filef = 'glg_tte_' + detector + '_'+f_list[0][2:]+f_list[1]+f_list[2]+'_'+f_list[3]
		file_list.append(filef)
	return file_list

def get_time_list(star,stop):
	t_star = int(time.Time(star).mjd*24)
	t_stop = time.Time(stop).mjd*24
	lis = np.arange(t_star,t_stop,1)/24
	time_list = time.Time(lis,format = 'mjd').fits
	return time_list
def open_daily_data_file(file_link):
	f = fits.open(file_link)
	time = f[2].data.field(0)
	ch = f[2].data.field(1)
	ch_n = f[1].data.field(0)
	e1 = f[1].data.field(1)
	e2 = f[1].data.field(2)
	return time,ch,ch_n,e1,e2

#time0 = ['2017-09-22 12:45:00','2017-09-22 22:00:00']
#t,ch,ch_n,e1,e2 = get_daily_data(time0,'n0')
#num = len(t)
#print(t[0],t[num-1])
#print(GBMtime.met_to_utc(t[0]).fits,GBMtime.met_to_utc(t[num-1]).fits)












