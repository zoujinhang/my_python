import numpy as np
#import astropy.units as u
from astropy.coordinates import cartesian_to_spherical,SkyCoord
from spherical_geometry.polygon import SphericalPolygon


class GBM_detector(object):
	def __init__(self,name,quaternion):
		self.name = name
		self.quaternion = quaternion

		self.position = self.get_position()
		#print(self.position)
		#self.gbm_xyz = np.array([0,0,1.0])
		p_lon = cartesian_to_spherical(self.position[0],self.position[1],self.position[2])[2].deg
		p_lat = cartesian_to_spherical(self.position[0],self.position[1],self.position[2])[1].deg
		#print(p_lon,p_lat)
		self.center = SkyCoord(ra = p_lon,dec = p_lat,frame = 'icrs',unit = 'deg')
		#print(self.center)

	def get_position(self):

		'''这里是计算探头在赤道坐标系中指向的过程。'''
		X = np.mat(self.gbm_xyz).T
		mat0 = self.get_mat(self.quaternion[0],self.quaternion[1],self.quaternion[2],self.quaternion[3])
		X1 = mat0*X
		x = np.array([X1[0],X1[1],X1[2]])
		x = np.array([x[0][0][0],x[1][0][0],x[2][0][0]])
		return x
	def get_mat(self,p1,p2,p3,p0):
		'''这里是计算出旋转矩阵'''
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
	def get_fov(self,radius):
		'''这里是计算出探头标注'''
		if radius >= 60:
			steps = 500
		elif radius >= 30:
			steps = 300
		else:
			steps = 250
		j2000 = self.center.icrs
		poly = SphericalPolygon.from_cone(j2000.ra.value,j2000.dec.value,radius,steps = steps)
		re =  [p for p in poly.to_radec()][0]
		return re
	def contains_point(self,point):

		steps = 300
		j2000 = self.center.icrs
		poly = SphericalPolygon.from_cone(j2000.ra.value,j2000.dec.value,self.radius,steps = steps)
		return poly.contains_point(point.cartesian.xyz.value)

class NaI0(GBM_detector):

	def __init__(self,quaternion,point = None):
		self.az = 45.89
		self.zen = 90 - 20.58
		self.radius = 60.0
		self.gbm_xyz = np.array([0.2446677589,0.2523893824,0.9361823057])
		super(NaI0, self).__init__('n0',quaternion)


class NaI1(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 45.11
		self.zen = 90 - 45.31
		self.radius = 60.0
		self.gbm_xyz = np.array([0.5017318971,0.5036621127,0.7032706462])
		super(NaI1, self).__init__('n1', quaternion)


class NaI2(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 58.44
		self.zen = 90 - 90.21
		self.radius = 60.0
		self.gbm_xyz = np.array([0.5233876659,0.8520868147,-0.0036651682])
		super(NaI2, self).__init__('n2', quaternion)


class NaI3(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 314.87
		self.zen = 90 - 45.24
		self.radius = 60.0
		self.gbm_xyz = np.array([0.5009495177,-0.5032279093,0.7041386753])
		super(NaI3, self).__init__('n3', quaternion)


class NaI4(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 303.15
		self.zen = 90. - 90.27
		self.radius = 60.0
		self.gbm_xyz = np.array([ 0.5468267487,-0.8372325378,-0.0047123847])
		super(NaI4, self).__init__('n4', quaternion)


class NaI5(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 3.35
		self.zen = 90 - 89.97
		self.radius = 60.0
		self.gbm_xyz = np.array([0.9982910766,0.0584352143,0.0005236008])
		super(NaI5, self).__init__('n5', quaternion)


class NaI6(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 224.93
		self.zen = 90 - 20.43
		self.radius = 60.0
		self.gbm_xyz = np.array([-0.2471260191,-0.2465229020,0.9370993606])
		super(NaI6, self).__init__('n6', quaternion)



class NaI7(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 224.62
		self.zen = 90 - 46.18
		self.radius = 60.0
		self.gbm_xyz = np.array([-0.5135631636,-0.5067957667,0.6923950822])
		super(NaI7, self).__init__('n7', quaternion)
class NaI8(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az =  236.61
		self.zen = 90 - 89.97
		self.radius = 60.0
		self.gbm_xyz = np.array([-0.5503349679,-0.8349438131,0.0005235846])
		super(NaI8, self).__init__('n8', quaternion)

class NaI9(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 135.19
		self.zen = 90 - 45.55
		self.radius = 60.0
		self.gbm_xyz = np.array([-0.5064476761,0.5030998708,0.7002865795])
		super(NaI9, self).__init__('n9', quaternion)


class NaIA(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 123.73
		self.zen = 90 - 90.42
		self.radius = 60.0
		self.gbm_xyz = np.array([-0.5552650628,0.8316411478,-0.0073303046])
		super(NaIA, self).__init__('na', quaternion)

class NaIB(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 183.74
		self.zen = 90 - 90.32
		self.radius = 60.0
		self.gbm_xyz = np.array([-0.9978547710,-0.0652279514,-0.0055850266])
		super(NaIB, self).__init__('nb', quaternion)

class BGO0(GBM_detector):

	def __init__(self,quaternion,point = None):
		self.az = 0.0
		self.zen = 0.0
		self.radius = 90.0
		self.gbm_xyz = np.array([1.0,0.0,0.0])
		super(BGO0, self).__init__('b0',quaternion)


class BGO1(GBM_detector):

	def __init__(self, quaternion, point=None):
		self.az = 180.0
		self.zen = 0.0
		self.radius = 90.0
		self.gbm_xyz = np.array([-1.0,0.0,0.0])
		super(BGO1, self).__init__('b1', quaternion)


#q = [0.09894184,0.81399423,0.56763536,0.07357984]
#n = NaI0(q)
#print(n.center)
