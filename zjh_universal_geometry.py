import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.table import Table
from astropy.coordinates import cartesian_to_spherical,SkyCoord,spherical_to_cartesian,get_sun,get_body_barycentric
from spherical_geometry.polygon import SphericalPolygon
import astropy.time as time
import astropy.units as u
from astropy.io import fits


class U_geometry(object):
	def __init__(self,quaternion,local_az,local_zen,sc_pos = None,time = None,ef_radius = 60):
		if time is not None:
			if (isinstance(time, str)):
				self.time = Time_transition.utc_time(time)
			else:
				self.time = Time_transition.met_to_utc(time)
		else:
			self.time = None
		self.radius = ef_radius
		self.quaternion = quaternion
		self.sc_pos = sc_pos
		self.local_az = local_az
		self.local_zen = local_zen
		self.detectors = Detectors(self.quaternion,self.local_az,self.local_zen)


	def get_detector_centers(self):
		return self.detectors.center_all
	def get_detector_index(self):
		return self.detectors.detector_index

	def get_good_detector_centers(self,source = None):
		good_detector_centers = []
		good_detector_index = []
		if source is not None:
			for index,center in enumerate(self.detectors.center_all):
				steps = 250
				j2000 = center.icrs
				poly = SphericalPolygon.from_cone(j2000.ra.value,j2000.dec.value,self.radius,steps = steps)
				if(poly.contains_point(source.cartesian.xyz.value)):
					good_detector_centers.append(center)
					good_detector_index.append(index)
		else:
			print('No source! return []')
		return good_detector_index,good_detector_centers
	def get_fov(self,radius):
		print('get_fov')
		return self.detectors.get_fov(radius)

	def get_separation(self,source = None):
		tab = Table(names=["Detector_index", "Separation"], dtype=["|S2", np.float64])
		if source is not None:
			for index,center in enumerate(self.detectors.center_all):
				sep = center.separation(source)
				tab.add_row([index,sep])
			tab['Separation'].unit = u.degree
			tab.sort('Separation')
			return tab
		else:
			print('No source! return None')
			return None
	def get_earth_point(self):
		if self.sc_pos is not None:
			self.calc_earth_points()
			return self.earth_points
		else:
			print('No satellite position!')
	def calc_earth_points(self):
		print('calc_earth_points')
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

	def detector_plot(self,radius = 10.0,source=None,points = None,good = False,projection = 'moll',lat_0 = 0,lon_0 = 180,ax = None,show_bodies = False):

		if ax is None:
			fig = plt.figure(figsize = (20,10))
			ax = fig.add_subplot(111)
		map = Basemap(projection=projection,lat_0=lat_0,lon_0 = lon_0,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax)

		fovs = self.get_fov(radius)
		if good and source :
			index,centor = self.get_good_detector_centers(source)
		else:
			index = self.get_detector_index()
			centor = self.get_detector_centers()
		if source:
			ra,dec = map(source.ra.value,source.dec.value)
			map.plot(ra,dec, '*', color='#f36c21' , markersize=20.)
		if points:
			ra,dec = map(points.ra.value,points.dec.value)
			map.plot(ra,dec, '*', color='#c7a252' , markersize=20.)
		for i in index:
			print(fovs[i])
			ra,dec = fovs[i]
			ra,dec = map(ra,dec)
			map.plot(ra,dec,'.',color = '#74787c',markersize = 3)
			x,y = map(centor[i].icrs.ra.value,centor[i].icrs.dec.value)
			plt.text(x-200000, y-200000,str(i), color='#74787c', size=22)
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
		if projection == 'moll':
			az1 = np.arange(0, 360, 30)
			zen1 = np.zeros(az1.size) + 2
			azname = []
			for i in az1:
				azname.append(r'${\/%s\/^{\circ}}$' % str(i))
			x1, y1 = map(az1, zen1)
			for index, value in enumerate(az1):
				plt.text(x1[index], y1[index], azname[index], size=20)
		map.drawmeridians(np.arange(0, 360, 30),dashes=[1,0],color='#d9d6c3')
		map.drawparallels(np.arange(-90, 90, 15),dashes=[1,0], labels=[1,0,0,1], color='#d9d6c3',size = 20)
		map.drawmapboundary(fill_color='#f6f5ec')



class Detectors(object):
	def __init__(self,quaternion,local_az,local_zen):
		self.quaternion = quaternion
		self.local_az = local_az
		self.local_zen = local_zen
		print(self.local_az,self.local_zen)
		self.vector = np.array(spherical_to_cartesian(1,self.local_az,self.local_zen)).T
		self.mat = self.get_mat(self.quaternion[0],self.quaternion[1],self.quaternion[2],self.quaternion[3])
		self.detector_index,self.center_all = self.get_center()

	def get_mat(self,p1,p2,p3,p0):
		mat = np.mat(np.zeros((3, 3)))
		mat[0, 0] = p0 ** 2 + p1 ** 2 - p2 ** 2 - p3 ** 2
		mat[0, 1] = 2 * (p1 * p2 - p0 * p3)
		mat[0, 2] = 2 * (p0 * p2 + p1 * p3)
		mat[1, 0] = 2 * (p3 * p0 + p2 * p1)
		mat[1, 1] = p0 ** 2 + p2 ** 2 - p3 ** 2 - p1 ** 2
		mat[1, 2] = 2 * (p2 * p3 - p1 * p0)
		mat[2, 0] = 2 * (p1 * p3 - p0 * p2)
		mat[2, 1] = 2 * (p0 * p1 + p3 * p2)
		mat[2, 2] = p0 ** 2 + p3 ** 2 - p1 ** 2 - p2 ** 2
		return mat

	def get_center(self):
		print('get_center')
		position_list = []
		detector_index_list = []

		for index,vector in enumerate(self.vector):
			X = np.mat(vector).T
			X1 = np.array(self.mat * X).T[0]
			position_list.append([X1[0],X1[1],X1[2]])
			detector_index_list.append(index)
		position_array = np.array(position_list).T
		detector_index = np.array(detector_index_list)
		position = cartesian_to_spherical(position_array[0],position_array[1],position_array[2])
		p_lon = position[2].deg
		p_lat = position[1].deg
		center_all = SkyCoord(ra = p_lon,dec = p_lat,frame = 'icrs',unit = 'deg')
		return detector_index,center_all

	def get_fov(self,radius):
		if radius >= 60:
			steps = 500
		elif(radius >=30):
			steps = 300
		else:
			steps = 250
		j2000 = self.center_all.icrs
		c = []
		for index,center in enumerate(j2000):
			print(str(index))
			poly = SphericalPolygon.from_cone(center.ra.value,center.dec.value,radius,steps=steps)
			re = [p for p in poly.to_radec()][0]
			#print(re)
			c.append(re)
		return c

class Time_transition(object):
	def __init__(self,start_time = None):
		self.start_time = 51910
		if(start_time is not None):
			print('Reminder : The start_time must be  astropy.time.Time() .')
			self.start_time = start_time.mjd

	@classmethod

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














