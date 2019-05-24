from scipy.stats.kde import gaussian_kde
from sklearn.neighbors import KernelDensity
from scipy.sparse import csc_matrix,eye,diags
from scipy.sparse.linalg import spsolve
import numpy as np

def density_map(tburst,photon,photon_max,photon_min,tburst_max,tburst_min, photon_bin, tburst_bin):
	Tresolution_sig= ((tburst_max - tburst_min)/tburst_bin)
	Eresolution_sig= ((photon_max - photon_min)/photon_bin)
	print('*************************************')
	print('this is burst infomation')
	print('time resolution is %s s' % Tresolution_sig)
	print('energy resolution is %s KeV' % Eresolution_sig)
	#利用高斯核估计密度分布
	index_sig = (photon >= photon_min) & (photon <= photon_max) & (tburst >= tburst_min) & (tburst <= tburst_max)

	tburst_sig=tburst[index_sig]
	photon_sig=photon[index_sig]
	total_sig=photon_sig.size
	print('total count in this region is %s' % total_sig)
	print('*************************************')
	k_sig = gaussian_kde(np.vstack([tburst_sig, photon_sig]))
	xi_sig, yi_sig = np.mgrid[tburst_min: tburst_max: Tresolution_sig,photon_min: photon_max: Eresolution_sig]
	density_sig = k_sig(np.vstack([xi_sig.flatten(), yi_sig.flatten()]))
	density_sig =density_sig.reshape(xi_sig.shape)
	return density_sig, xi_sig, yi_sig, total_sig
#baseline_utils:
def WhittakerSmooth(x,w,lambda_):
	X=np.matrix(x)#这里将数组转化为矩阵。矩阵之后就不可以用索引进行引用了。
	m=X.size
	i=np.arange(0,m)
	E=eye(m,format='csc')
	D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
	W=diags(w,0,shape=(m,m))
	A=csc_matrix(W+(lambda_*D.T*D))
	B=csc_matrix(W*X.T)
	background=spsolve(A,B)
	return np.array(background)
def airPLS(x, lambda_=100, itermax=15):
	m=x.shape[0]
	w=np.ones(m)
	for i in range(1,itermax+1):
		z=WhittakerSmooth(x,w,lambda_)
		d=x-z
		dssn=np.abs(d[d<0].sum())
		if(dssn<0.001*(abs(x)).sum() or i==itermax):
			if(i==itermax):
				print('WARING max iteration reached!')
			break
		w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
		w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
		w[0]=np.exp(i*(d[d<0]).max()/dssn)
		w[-1]=w[0]
	return z
def get_contour_verts(cn, level):
	contours = []
	cc=cn.collections[level]
	paths = []
	for pp in cc.get_paths():
		xy = []
		for vv in pp.iter_segments():
			xy.append(vv[0])
		paths.append(np.vstack(xy))
	contours.append(paths)
	contours = np.vstack(contours)
	return contours
'''
def get_contour_verts_new(cn):
	contours = []
	# for each contour line
	for cc in cn.collections:
		paths = []
		# for each separate section of the contour line
		for pp in cc.get_paths():
			xy = []
			# for each segment of that section
        		for vv in pp.iter_segments():
				xy.append(vv[0])
			paths.append(np.vstack(xy))
		contours.append(paths)

	return contours

'''






