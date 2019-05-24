from zjh_gbm_detector import *
from zjh_gbm_time import *
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import astropy.units as u
from astropy.table import Table
from mpl_toolkits.basemap import Basemap
from astropy.coordinates import get_sun,get_body_barycentric,cartesian_to_spherical
from astropy.coordinates import SkyCoord

class GBM(object):
	'''
	这里是GBM的主函数
	'''
	def __init__(self,quaternion,sc_pos = None,gbm_time = None):
		if gbm_time is not None:
			if (isinstance(gbm_time, str)):
				self.time = GBMtime.utc_time(gbm_time)
			else:
				self.time = GBMtime.met_to_utc(gbm_time)
		else:
			self.time = None
		self.n0 = NaI0(quaternion)
		self.n1 = NaI1(quaternion)
		self.n2 = NaI2(quaternion)
		self.n3 = NaI3(quaternion)
		self.n4 = NaI4(quaternion)
		self.n5 = NaI5(quaternion)
		self.n6 = NaI6(quaternion)
		self.n7 = NaI7(quaternion)
		self.n8 = NaI8(quaternion)
		self.n9 = NaI9(quaternion)
		self.na = NaIA(quaternion)
		self.nb = NaIB(quaternion)
		self.b0 = BGO0(quaternion)
		self.b1 = BGO1(quaternion)
		self.detectors = OrderedDict(n0=self.n0,n1=self.n1,n2=self.n2,n3=self.n3,n4=self.n4,n5=self.n5,n6=self.n6,n7=self.n7,n8=self.n8,n9=self.n9,na=self.na,nb=self.nb,b0=self.b0,b1=self.b1)
		self.NaI_detectors = OrderedDict(n0=self.n0,n1=self.n1,n2=self.n2,n3=self.n3,n4=self.n4,n5=self.n5,n6=self.n6,n7=self.n7,n8=self.n8,n9=self.n9,na=self.na,nb=self.nb)
		self.BGO_detectors = OrderedDict(b0=self.b0,b1=self.b1)
		self.quaternion = quaternion
		self.sc_pos = sc_pos


	def get_centers(self,NaI = True,BGO = True):
		'''
		这里是获取探头指向中心的过程函数可以通过以下参数设置需要的探头种类。
		:param NaI: 选择 NaI 探头，默认为选择。
		:param BGO: 选择 BGO 探头，默认为选择。
		:return: 返回以探头名称为索引的字典。
		'''

		centers = {}
		if(NaI):
			for key in self.NaI_detectors.keys():
				centers[self.NaI_detectors[key].name] = self.NaI_detectors[key].center
		if(BGO):
			for key in self.BGO_detectors.keys():
				centers[self.BGO_detectors[key].name] = self.BGO_detectors[key].center
		return centers

	def get_good_centers(self,point = None,NaI = True,BGO = True):
		'''
		获取有效探头的指向中心。
		:param point: 源的指向点。
		:param NaI: 选择探头种类。
		:param BGO: 选择探头种类。
		:return: 返回以探头名称为索引的字典。
		'''

		if point is not None:
			centers = {}
			if(NaI):
				for key in self.NaI_detectors.keys():
					if(self.NaI_detectors[key].contaions_point(point)):
						centers[self.NaI_detectors[key].name] = self.NaI_detectors[key].center
			if(BGO):
				for key in self.BGO_detectors.keys():
					if(self.BGO_detectors[key].contaions_point(point)):
						centers[self.BGO_detectors[key].name] = self.BGO_detectors[key].center
			return centers
		else:
			print('没有源，无法确定有效探头！')
			return False

	def get_fov(self,radius,NaI = True,BGO = True):

		polys = {}
		detector_list = []
		if(NaI):
			for key in self.NaI_detectors.keys():
				polys[self.NaI_detectors[key].name] = self.NaI_detectors[key].get_fov(radius)
				detector_list.append(key)
		if(BGO):
			for key in self.BGO_detectors.keys():
				polys[self.BGO_detectors[key].name] = self.BGO_detectors[key].get_fov(radius)
				detector_list.append(key)
		return detector_list,polys

	def get_good_fov(self,radius,point = None,NaI = True,BGO = True):
		if point is not None:
			polys = {}
			detector_list = []
			if(NaI):
				for key in self.NaI_detectors.keys():
					if(self.NaI_detectors[key].contains_point(point)):
						polys[self.NaI_detectors[key].name] = self.NaI_detectors[key].get_fov(radius)
						detector_list.append(key)
			if(BGO):
				for key in self.BGO_detectors.keys():
					if(self.BGO_detectors[key].contains_point(point)):
						polys[self.BGO_detectors[key].name] = self.BGO_detectors[key].get_fov(radius)
						detector_list.append(key)

			return detector_list,polys
		else:
			print('没有源，无法确定有效探头！')
			return False

	def get_separation(self,source = None,NaI = True,BGO = True):

		tab = Table(names=["Detector", "Separation"], dtype=["|S2", np.float64])
		if source is not None:
			if(NaI):
				for key in self.NaI_detectors.keys():
					sep = self.NaI_detectors[key].center.separation(source)
					tab.add_row([key, sep])
			if(BGO):
				for key in self.BGO_detectors.keys():
					sep = self.BGO_detectors[key].center.separation(source)
					tab.add_row([key, sep])
			tab['Separation'].unit = u.degree
			tab.sort("Separation")
			return tab
		else:
			print('没有源，无法计算夹角！')
	def get_earth_point(self):
		if self.sc_pos is not None:
			self.calc_earth_points()
			return self.earth_points
		else:
			print('没有卫星位置！')

	def calc_earth_points(self):
		'''
		计算地球遮挡区域，无返回值。
		:return:
		'''
		xyz_position = SkyCoord(x=self.sc_pos[0],y=self.sc_pos[1],z=self.sc_pos[2],frame='icrs',representation='cartesian')
		earth_radius = 6371. * u.km
		fermi_radius = np.sqrt((self.sc_pos ** 2).sum())
		horizon_angle = 90 - np.rad2deg(np.arccos((earth_radius / fermi_radius).to(u.dimensionless_unscaled)).value)
		horizon_angle = (180 - horizon_angle) * u.degree
		num_points = 3000
		ra_grid_tmp = np.linspace(0, 360, num_points)
		dec_range = [-90, 90]
		cosdec_min = np.cos(np.deg2rad(90.0 + dec_range[0]))
		cosdec_max = np.cos(np.deg2rad(90.0 + dec_range[1]))
		v = np.linspace(cosdec_min, cosdec_max, num_points)
		v = np.arccos(v)
		v = np.rad2deg(v)
		v -= 90.
		dec_grid_tmp = v
		ra_grid = np.zeros(num_points ** 2)
		dec_grid = np.zeros(num_points ** 2)
		itr = 0
		for ra in ra_grid_tmp:
			for dec in dec_grid_tmp:
				ra_grid[itr] = ra
				dec_grid[itr] = dec
				itr += 1
		all_sky = SkyCoord(ra=ra_grid, dec=dec_grid, frame='icrs', unit='deg')
		condition = all_sky.separation(xyz_position) > horizon_angle
		self.earth_points = all_sky[condition]




	def detector_plot(self,radius = 60.0,point = None,good = False,projection = 'moll',lat_0 = 0,lon_0 = 180,map = None, NaI = True,BGO = True,show_bodies = False):
		
		map_flag = False

		if map is None:
			fig = plt.figure(figsize = (20,20))
			ax = fig.add_subplot(111)
			map = Basemap(projection=projection,lat_0=lat_0,lon_0 = lon_0,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax)
		else:
			map_flag = True

		if good and point:
			detector_list,fovs = self.get_good_fov(radius = radius,point = point,NaI = NaI,BGO = BGO)

		else:
			detector_list,fovs = self.get_fov(radius,NaI = NaI,BGO = BGO)
		if point:
			ra,dec = map(point.ra.value,point.dec.value)
			map.plot(ra,dec , '*', color='#f36c21' , markersize=20.)

		for key in detector_list:
			ra,dec = fovs[self.detectors[key].name]
			ra,dec = map(ra,dec)
			map.plot(ra,dec,'.',color = '#74787c',markersize = 3)
			x,y = map(self.detectors[key].center.icrs.ra.value,self.detectors[key].center.icrs.dec.value)
			plt.text(x-200000, y-200000,self.detectors[key].name, color='#74787c', size=22)
		if show_bodies and self.sc_pos is not None:
			earth_points = self.get_earth_point()
			lon, lat = earth_points.ra.value, earth_points.dec.value
			lon,lat = map(lon,lat)
			map.plot(lon, lat, ',', color="#0C81F9", alpha=0.1, markersize=4.5)
			if self.time is not None:
				earth_r = get_body_barycentric('earth',self.time)
				moon_r = get_body_barycentric('moon',self.time)
				r_e_m = moon_r - earth_r
				r = self.sc_pos -np.array([r_e_m.x.value,r_e_m.y.value,r_e_m.z.value])*u.km
				moon_point_d = cartesian_to_spherical(-r[0],-r[1],-r[2])
				moon_ra,moon_dec = moon_point_d[2].deg,moon_point_d[1].deg
				moon_point = SkyCoord(moon_ra,moon_dec,frame='icrs', unit='deg')
				moon_ra,moon_dec = map(moon_point.ra.deg,moon_point.dec.deg)
				map.plot(moon_ra,moon_dec,'o',color = '#72777b',markersize = 20)
				plt.text(moon_ra,moon_dec-800000,'moon',size = 20)
			if show_bodies and self.time is not None:
				tmp_sun = get_sun(self.time)
				sun_position = SkyCoord(tmp_sun.ra.deg,tmp_sun.dec.deg,unit='deg', frame='icrs')
				sun_ra,sun_dec = map(sun_position.ra.value,sun_position.dec.value)
				map.plot(sun_ra,sun_dec ,'o',color = '#ffd400',markersize = 40)
				plt.text(sun_ra-550000,sun_dec-200000,'sun',size = 20)


		if not map_flag:
			if projection == 'moll':
				az1 = np.arange(0, 360, 30)
				zen1 = np.zeros(az1.size) + 2
				azname = []
				for i in az1:
					azname.append(r'${\/%s\/^{\circ}}$' % str(i))
				x1, y1 = map(az1, zen1)
				for index, value in enumerate(az1):
					plt.text(x1[index], y1[index], azname[index], size=20)
			_ = map.drawmeridians(np.arange(0, 360, 30),dashes=[1,0],color='#d9d6c3')
			_ = map.drawparallels(np.arange(-90, 90, 15),dashes=[1,0], labels=[1,0,0,1], color='#d9d6c3',size = 20)
			map.drawmapboundary(fill_color='#f6f5ec')



#q = [0.09894184,0.81399423,0.56763536,0.07357984]
#myGBM = GBM(q)
#tb = myGBM.NaI_detectors.keys()
#print(list(tb))




