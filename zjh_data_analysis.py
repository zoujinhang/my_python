#import zzh_py3_r_baseline as zhbl
import rpy2.robjects as robjects
from rpy2.robjects import r
#from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import numpy as np
from astropy.stats import bayesian_blocks
import os
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit

robjects.r("library(baseline)")
robjects.numpy2ri.activate()
r("source('/home/laojin/software/my_r/zjh_r_remv_photo.R')")

def get_all_energy(time,ch,ch_n,e1,e2):
	new_t = np.array([])
	new_energy = np.array([])
	for index,channel in enumerate(ch_n):
		ch_t_index = np.where(ch == channel)
		ch_t = time[ch_t_index]
		energy_array = get_energy_of_ch(ch_t,e1[index],e2[index])
		new_t = np.concatenate((new_t,ch_t))
		new_energy = np.concatenate((new_energy,energy_array))
	index_all = np.argsort(new_t)
	new_t = new_t[index_all]
	new_energy = new_energy[index_all]
	return new_t,new_energy

def r_remv_photo(time,ch,ch_n,bin_size=1.,bins = None,inti = None,hwi = None):
	new_t = np.array([])
	new_ch = np.array([])
	if(bins == None):
		t_bin_array = np.arange(time[0],time[-1]+bin_size,bin_size)
	else:
		t_bin_array = bins
	if(inti == None):
		inti = int(len(t_bin_array)/10)
	if(hwi == None):
		hwi = int(20/bin_size)
	for index,value in enumerate(ch_n):
		ch_t_index = np.where(ch == value)
		ch_t = time[ch_t_index]
		bin_n,bin_edges = np.histogram(ch_t,bins = t_bin_array)
		ch_tc = (bin_edges[1:]+bin_edges[:-1])/2
		ch_tc,bin_n,back = r_baseline(ch_tc,bin_n,hwi=hwi,inti = inti)
		ch_tc, bin_n1, back = r_baseline(ch_tc,np.abs(bin_n))#,lam = 2,hwi = 500,it = 10
		bin_n = np.abs(bin_n) - back.sum()/back.size
		index3 = np.where(bin_n < 0)
		if(len(bin_n[index3])>0):
			nnn = bin_n[index3].sum() / (len(bin_n[index3]))
			bin_n = np.abs(bin_n) + 1.5*nnn
		else:
			bin_n = np.abs(bin_n)
		#bin_n = np.abs(bin_n)-back
		ch_tc, bin_n1, back = r_baseline(ch_tc,np.abs(bin_n))
		#ch_tc, bin_n, back = r_baseline(ch_tc, np.abs(bin_n), lam=3, hwi=300, it=7)
		#bin_n = np.abs(bin_n)-back
		bin_n = np.abs(bin_n) - back.sum()/back.size
		index3 = np.where(bin_n < 0)
		if(len(bin_n[index3])>0):
			nnn = bin_n[index3].sum() / (len(bin_n[index3]))
			bin_n = np.abs(bin_n) + 1.5*nnn
		else:
			bin_n = np.abs(bin_n)
		n_ch_t = np.array([])
		for j,i in enumerate(ch_tc):
			b_t_index = np.where((ch_t >= bin_edges[j]) & (ch_t < (bin_edges[j+1])))
			b_t = ch_t[b_t_index]

			rand_array = np.random.rand(len(b_t))
			rand_index = np.argsort(rand_array)
			eee = int(np.abs(bin_n[j]))

			if (eee > 0):
				index_c = rand_index[0:eee]
				b_t_c = b_t[index_c]
				n_ch_t = np.concatenate((n_ch_t,b_t_c))
		n_ch_t_index = np.argsort(n_ch_t)
		n_ch_t = n_ch_t[n_ch_t_index]
		b_ch = np.full(len(n_ch_t),value)
		new_t = np.concatenate((new_t,n_ch_t))
		new_ch = np.concatenate((new_ch,b_ch))
	new_index = np.argsort(new_t)
	new_t = new_t[new_index]
	new_ch = new_ch[new_index]
	return new_t,new_ch

def fast_remv_photo(time,ch,ch_n,bin_size = 1.,bins = None,inti = None,hwi = None):
	if(bins == None):
		t_bin_array = np.arange(time[0],time[-1]+bin_size,bin_size)
	else:
		t_bin_array = bins
	if(inti == None):
		inti = int(len(t_bin_array)/10)
	if(hwi == None):
		hwi = int(60/bin_size)
	r.assign('ch', ch*1.)
	r.assign('t',time)
	r.assign('ch_n',ch_n)
	r.assign('bin_size',bin_size)
	r.assign('bins',t_bin_array)
	r.assign('inti',inti)
	r.assign('hwi',hwi)
	r("r_data = remv_photo(t,ch,ch_n,bins,bin_size=bin_size ,inti=inti,hwi=hwi)")
	r('new_t = r_data$new_t')
	r('new_ch = r_data$new_ch')
	new_t = np.float32(r('new_t'))
	new_ch = np.float32(r('new_ch'))
	return new_t,new_ch

def get_energy_of_ch(time,e1,e2):
	numb = len(time)
	energy_random_arr = np.random.random_sample(numb)
	energy_array = e1 + (e2-e1)*energy_random_arr
	return energy_array

def r_baseline(time,rate,lam = None,hwi = None,it = None,inti = None):
	r.assign('rrate',rate)
	r('y = matrix(rrate,nrow = 1)')
	dt = time[1]-time[0]
	#fillpeak_int = str(int(time[-1]-time[0])/dt/10)
	if(lam is None):
		lam = int(0.708711*(len(rate))**0.28228+0.27114)
	if(hwi is None):
		hwi = int(40/dt)
	if(it is None):
		it = 10

	if(inti is None):
		fillpeak_int = str(int(len(rate)/10))
	else:
		fillpeak_int = str(inti)

	if(lam < 1):
		lam = 1
	r("rbase = baseline(y,lam ="+str(lam)+",hwi = "+str(hwi)+",it = "+str(it)+",int = "+fillpeak_int+",method = 'fillPeaks')")
	r("bs = getBaseline(rbase)")
	r("cs = getCorrected(rbase)")
	bs = np.float32(r('bs')[0])
	cs = np.float32(r('cs')[0])
	return time,cs,bs

def easy_histogram(time,binsize = 1.):

	tbin_array = np.arange(time[0],time[-1]+binsize,binsize)
	bin_n,bin_edges = np.histogram(time,bins=tbin_array)
	x = (bin_edges[1:]+bin_edges[:-1])/2
	return x,bin_n/binsize


def get_peak(t,value,order = 12,num = 10):

	max_index = np.argmax(value)
	t_max = t[max_index]
	index = np.where((t>=t_max-7)&(t<=t_max+7))
	t_p = t[index]
	value_p = value[index]
	fit = np.polyfit(t_p, value_p, order)
	t_n = np.linspace(t_p[0], t_p[-1], 1000)
	v_n = np.polyval(fit,t_n)

	max_index = np.argmax(v_n)
	t_max = t_n[max_index]

	t_array = np.zeros(num)
	for i in range(num):
		index = np.where((t>=t_max-7)&(t<=t_max+7))
		t_p = t[index]
		value_p = value[index]
		fit = np.polyfit(t_p, value_p, order)
		t_n = np.linspace(t_p[0], t_p[-1], 1000)
		v_n = np.polyval(fit,t_n)
		max_index = np.argmax(v_n)
		t_max = t_n[max_index]
		t_array[i] = t_max

	t_peak = t_array.sum()/t_array.size
	return t_peak

def get_lag(x1,x2,dt,x1_err = None,x2_err = None,num = 1000,order = 12,name = None,savedir = None):

	if(x1_err is None):
		x1_err = np.sqrt(x1)
	if(x2_err is None):
		x2_err = np.sqrt(x2)
	n = len(x1)
	lag_t = np.arange(-n+1,n,1)*dt
	nccf = np.correlate(x1/x1.max(),x2/x2.max(),'full')
	lag_ever = get_peak(lag_t,nccf,order=order)
	print('--------------------------------------------------')

	lag_list = []
	if((name is not None)and(savedir is not None)):
		print('We plot the check picture in '+savedir)
		print('The name is ' + name+'.png')
		index_array = np.where((lag_t >= lag_ever-7)&(lag_t<=lag_ever+7))
		fit_nccf = nccf[index_array]
		fit_lag_t = lag_t[index_array]
		new_t = np.linspace(fit_lag_t[0],fit_lag_t[-1],1000)
		c = np.polyfit(fit_lag_t,fit_nccf,order)
		yy = np.polyval(c,new_t)
		if(os.path.exists(savedir) == False):
			os.makedirs(savedir)
		plt.plot(fit_lag_t,fit_nccf)
		plt.plot(new_t,yy)
		plt.savefig(savedir+name+'.png')
		plt.close()

	for i in range(num):
		rand = np.random.randn(n)
		x1_sim = x1 + x1_err*rand
		x2_sim = x2 + x2_err*rand
		nccf = np.correlate(x1_sim/x1_sim.max(),x2_sim/x2_sim.max(),'full')
		lag = get_peak(lag_t,nccf,order = order,num = 5)
		lag_list.append(lag)
		print('Monte Carlo Lag '+str(i)+' :',lag,end='\r')

	lag_array = np.array(lag_list)
	lag_err = lag_array.std()
	print('The measured lag:', lag_ever)
	print('The lag error', lag_err)
	print('--------------------------------------------------')
	return lag_ever,lag_err


def check_rate(rate,standard = 10):
	check = 0
	n = len(rate)-1

	for index in range(1,n):
		if((rate[index] == 0)&(rate[index-1] == 0)&(rate[index+1] == 0)):
			check = check+1
			if(check > standard):
				return False
		else:
			check = 0
	return True

def get_pulse_duration(t,edges = None,cafe_time = 0,r = 3,savedir = None,lam = None,hwi = int(20/0.02),gamma = None,p0 = 0.005):
	'''
	:param t: data time of point
	:param edges: [time_start,time_stop]
	:param cafe_time:
	:param savedir:
	:return: start_edges,stop_edges
	'''
	if(edges is not None):
		bins = np.arange(edges[0],edges[-1]+0.02,0.02)
		bin_n,bin_edges = np.histogram(t,bins=bins)
		t_h = (bin_edges[1:]+bin_edges[:-1])/2
		rate = bin_n
	else:
		t_h,rate = easy_histogram(t,binsize = 0.02)
		index = np.where((t_h>t_h[0]+5)&(t_h<t_h[-1]-5))
		t_h = t_h[index]
		rate = rate[index]*0.02
	start_edges = []
	stop_edges = []
	if(check_rate(rate) == False):
		print('the data is not good we can not get the pulse`s duration.')
		return start_edges,stop_edges

	print('getting pulse duration...')
	t_r,cs,bs = r_baseline(t_h,rate,lam = lam,hwi = hwi)
	print('r_baseline useing parameter are default.')
	bs_ev = bs.sum()/bs.size
	va = cs + bs_ev
	edges = bayesian_blocks(t_r,np.round(va),fitness = 'events',gamma = gamma,p0 = p0)
	#edges = bayesian_blocks(t_r,va,np.sqrt(bs_ev),fitness = 'measures',p0 = p0,gamma = gamma)
	if(len(edges)>3):
		print('step 1 ....')
		if(savedir is not None):
			if(os.path.exists(savedir) == False):
				os.makedirs(savedir)
			plt.plot(t_h,rate)
			plt.plot(t_r,bs)
			plt.savefig(savedir+'check_1_step.png')
			plt.close()
		binsize = edges[1:]-edges[:-1]
		binstart = edges[:-1]
		binstop = edges[1:]
		bin_high = []
		#bin_high_sigma = []
		for index1, value1 in enumerate(binstart):
			index_in_bin = np.where((t_h >= value1) & (t_h < binstop[index1]))
			bs_bin = va[index_in_bin]
			#bin_high_sigma.append(np.std(bs_bin))
			bin_high.append(np.mean(bs_bin))
		bin_high = np.array(bin_high)
		#bin_high_sigma = np.array(bin_high_sigma)
		trait = []
		for index1, hight in enumerate(bin_high):
			if (index1 == 0):
				trait.append('start')
			elif (index1 == len(bin_high) - 1):
				trait.append('stop')
			else:
				if ((hight > bin_high[index1 - 1]) and (hight > bin_high[index1 + 1])):
					trait.append('pulse')
				elif ((hight < bin_high[index1 - 1]) and (hight < bin_high[index1 + 1])):
					trait.append('cafe')
				else:
					trait.append('fringe')
		sort_index = np.argsort(binsize)
		sort_binsize = binsize[sort_index]
		sort_bin_hight = bin_high[sort_index]
		#sort_bin_hight_sigma = bin_high_sigma[sort_index]
		bin_big3 = sort_binsize[-3:]
		hight_big3 = sort_bin_hight[-3:]
		#hight_sigma_big3 = sort_bin_hight_sigma[-3:]
		eva_hight = (bin_big3*hight_big3).sum()/bin_big3.sum()
		#eva_hight_sigma = (hight_sigma_big3*bin_big3).sum()/bin_big3.sum()
		print('step 2 ....')
		if (savedir is not None):
			binhigh = np.concatenate(([bin_high[0]], bin_high))
			#binhigh_sigma = np.concatenate(([bin_high_sigma[0]], bin_high_sigma))
			plt.plot(t_h, va)
			for edge in edges:
				plt.axvline(x=edge, color='r')
			#plt.plot(edges,binhigh+binhigh_sigma,linestyle = 'steps',color ='g')
			#plt.plot(edges,binhigh-binhigh_sigma,linestyle = 'steps',color ='g')
			plt.plot(edges, binhigh, linestyle='steps', color='k')
			plt.plot([t_h[0],t_h[-1]],[eva_hight,eva_hight],linestyle = 'steps',color ='y')
			plt.plot([t_h[0],t_h[-1]],[eva_hight-r,eva_hight-r],linestyle = 'steps',color ='y')
			plt.plot([t_h[0],t_h[-1]],[eva_hight+r,eva_hight+r],linestyle = 'steps',color ='y')
			#plt.plot([t_h[0],t_h[-1]],[eva_hight+3*eva_hight_sigma,eva_hight+r],linestyle = 'steps',color ='y')

			plt.savefig(savedir + 'check_2_step.png')
			plt.close()


		binstart1 = [binstart[0]]
		binstop1 = []
		for index1,tr in enumerate(trait):
			if((tr != 'start')and(tr != 'stop')):
				if((bin_high[index1-1]>eva_hight-r)and(bin_high[index1+1]>eva_hight - r)and(bin_high[index1]>eva_hight - r)):
					binstart1.append(binstart[index1])
					binstop1.append(binstop[index1])
		binstop1.append(binstop[-1])
		edges1 = np.array(binstart1+binstop1)
		edges1 = np.unique(edges1)
		edges1 = np.sort(edges1)
		binstart1 = edges1[:-1]
		binstop1 = edges1[1:]
		binsize1 = binstop1 - binstart1
		bin_high1 = []
		#bin_high_sigma1 = []
		for index1, value1 in enumerate(binstart1):
			index_in_bin = np.where((t_h >= value1) & (t_h < binstop1[index1]))
			bs_bin = va[index_in_bin]
			#bin_high_sigma1.append(np.std(bs_bin))
			bin_high1.append(np.mean(bs_bin))
		bin_high1 = np.array(bin_high1)
		#bin_high_sigma1 = np.array(bin_high_sigma1)
		trait1 = []
		for index1,hight in enumerate(bin_high1):
			if(index1 == 0):
				trait1.append('start')
			elif(index1 == len(bin_high1)-1):
				trait1.append('stop')
			else:
				if((hight>bin_high1[index1-1])and(hight>bin_high1[index1+1])):
					trait1.append('pulse')
				elif((hight<bin_high1[index1-1])and(hight<bin_high1[index1+1])):
					trait1.append('cafe')
				else:
					trait1.append('fringe')
		print('step 3 ....')
		if (savedir is not None):
			binhigh1 = np.concatenate(([bin_high1[0]], bin_high1))
			#binhigh_sigma1 = np.concatenate(([bin_high_sigma1[0]], bin_high_sigma1))
			plt.plot(t_h,va)
			for edge in edges1:
				plt.axvline(x=edge, color='r')
			#plt.plot(edges1, binhigh1-binhigh_sigma1, linestyle='steps', color='g')
			#plt.plot(edges1, binhigh1+binhigh_sigma1, linestyle='steps', color='g')
			plt.plot(edges1, binhigh1, linestyle='steps', color='k')
			plt.savefig(savedir + 'check_3_step.png')
			plt.close()

		#------------------------------------------------------

		#--------------------------------------------------------
		binstart3 = [binstart1[0]]
		binstop3 = []
		for index1,tr in enumerate(trait1):
			if((binsize1[index1]<=30)and(bin_high1[index1]>=eva_hight+r)):
				binstart3.append(binstart1[index1])
				binstop3.append(binstop1[index1])
		binstop3.append(binstop1[-1])
		edges3 = np.array(binstart3+binstop3)
		edges3 = np.unique(edges3)
		edges3 = np.sort(edges3)
		binstart3 = edges3[:-1]
		binstop3 = edges3[1:]
		binsize3 = binstop3 - binstart3
		bin_high3 = []
		#bin_high_sigma3 = []
		for index1, value1 in enumerate(binstart3):
			index_in_bin = np.where((t_h >= value1) & (t_h < binstop3[index1]))
			bs_bin = va[index_in_bin]
			#bin_high_sigma3.append(np.std(bs_bin))
			bin_high3.append(np.mean(bs_bin))
		bin_high3 = np.array(bin_high3)
		#bin_high_sigma3 = np.array(bin_high_sigma3)
		trait3 = []
		for index1,hight in enumerate(bin_high3):
			if(index1 == 0):
				trait3.append('start')
			elif(index1 == len(bin_high3)-1):
				trait3.append('stop')
			else:
				if((hight>bin_high3[index1-1])and(hight>bin_high3[index1+1])):
					trait3.append('pulse')
				elif((hight<bin_high3[index1-1])and(hight<bin_high3[index1+1])):
					trait3.append('cafe')
				else:
					trait3.append('fringe')
		#------------------------------------------------------

		print('step 4 ....')
		if (savedir is not None):
			binhigh3 = np.concatenate(([bin_high3[0]], bin_high3))
			#binhigh_sigma3 = np.concatenate(([bin_high_sigma3[0]], bin_high_sigma3))
			plt.plot(t_h,va)
			for edge in edges3:
				plt.axvline(x=edge, color='r')
			#plt.plot(edges3, binhigh3-binhigh_sigma3, linestyle='steps', color='g')
			#plt.plot(edges3, binhigh3+binhigh_sigma3, linestyle='steps', color='g')
			plt.plot(edges3, binhigh3, linestyle='steps', color='k')
			plt.savefig(savedir + 'check_4_step.png')
			plt.close()
		#--------------------------------------------------------
		start_edges4 = []
		stop_edges4 = []
		print('step 5 ....')
		if(savedir is not None):
			plt.plot(t_h,rate)
		if(len(edges3)>3):
			start_edges4.append(binstart3[1])
			for i in range(1,len(binstart3)-1):
				size = binsize3[i]
				if ((size > cafe_time)and(i != len(binstart3)-2)):
					if((bin_high3[i-1] > bin_high3[i]) & (bin_high3[i] < bin_high3[i+1])):
						start_edges4.append(binstop3[i])
						stop_edges4.append(binstart3[i])
			stop_edges4.append(binstop3[-2])
			start_edges4 = np.array(start_edges4)
			stop_edges4 = np.array(stop_edges4)
			start_edges_c =[]
			stop_edges_c = []
			for i in range(len(start_edges4)):
				if(start_edges4[i] != stop_edges4[i]):
					start_edges_c.append(start_edges4[i])
					stop_edges_c.append(stop_edges4[i])
			start_edges = np.array(start_edges_c)
			stop_edges = np.array(stop_edges_c)
			if(savedir is not None):
				for index1,value1 in enumerate(start_edges):
					plt.axvline(x = value1,color = 'r')
					plt.axvline(x = stop_edges[index1],color = 'g')
		if(savedir is not None):
			plt.savefig(savedir+'check_5_step.png')
			plt.close()
	return start_edges,stop_edges


def get_sigma(x,p,standard = 0.94):
	'''

	:param x:
	:param p: 一定是概率,p.sum() = 1
	:param standard:
	:return:
	'''
	n = len(x)
	k = 0
	xk = []
	x = np.array(x)
	p = np.array(p)
	indexmax = np.argmax(p)
	index_r = indexmax+1
	index_l = indexmax-1
	pk = p[indexmax]
	xk.append(x[indexmax])
	while (pk<standard)and(k<n):
		if((index_r<n-1)and(index_l>0)):
			if(p[index_r]<p[index_l]):
				xk.append(x[index_l])
				pk = pk + p[index_l]
				index_l = index_l - 1
			else:

				xk.append(x[index_r])
				pk = pk + p[index_r]
				index_r = index_r + 1

		elif((index_r >= n-1)and(index_l>0)):
			xk.append(x[index_l])
			pk = pk + p[index_l]
			index_l = index_l - 1
		elif((index_l <= 0)and(index_r<n)):
			xk.append(x[index_r])
			pk = pk + p[index_r]
			index_r = index_r + 1
		k = k+1
	x_t = np.max(xk)
	x_b = np.min(xk)
	length = (x_t-x_b)*0.5
	local = (x_t+x_b)*0.5
	return local,length

def get_params(x,y,E0,z,sigma = None,n=1):
	R_M = 0.315
	R_A = 1-R_M
	value = quad(lambda x:(1+x)**n/np.sqrt(R_M*(1+x)**3+R_A),0,z)[0]
	def function1(e, t, a, eq):
		return t * (e ** a - E0 ** a) - (1 + n) / 2 / 67.3 * (e ** n - E0 ** n) / eq**n * value
	p,l = curve_fit(function1,x,y,method='dogbox',bounds=([0,0,1e-10], [np.inf,np.inf,np.inf]),sigma=sigma)
	return p,function1

def fit_lags_isotropic(e,lag,lag_err,n,E0,z,num = 5000):
	e = np.array(e)
	lag = np.array(lag)
	lag_err = np.array(lag_err)
	p,function1 = get_params(e,lag,E0,z,sigma=lag_err,n = n)
	t = []
	a = []
	eq = []
	while num:
		try:
			lag_sim = lag + lag_err*np.random.randn(lag_err.size)
			p_sim,function2 = get_params(e,lag_sim,E0,z,sigma = lag_err,n = n)
			t.append(p_sim[0])
			a.append(p_sim[1])
			eq.append(p_sim[2])
			num = num-1
			print('+++',end = '\r')
		except (RuntimeError):
			print('---',end = '\r')

	t = np.array(t)
	a = np.array(a)
	eq = np.array(eq)
	logeq = np.log10(eq)
	print('....')
	t_edges = np.arange(t.min(),t.max()+0.01,0.01)
	bin_n_t,bin_edge_t = np.histogram(t,bins=t_edges)
	bin_edge_t_c = (bin_edge_t[1:]+bin_edge_t[:-1])*0.5
	t_local,t_err_len = get_sigma(bin_edge_t_c,bin_n_t/bin_n_t.sum(),standard=0.94)
	t_start = t_local-t_err_len
	t_stop = t_local + t_err_len
	t_local1,t_err_len1 = get_sigma(bin_edge_t_c,bin_n_t/bin_n_t.sum(),standard=0.68)
	t_err_l = p[0]-t_local1+t_err_len1
	t_err_h = t_local1+t_err_len1-p[0]
	t_edges = np.linspace(t_start,t_stop,50)
	bin_n_t,bin_edge_t = np.histogram(t,bins=t_edges)
	bin_n_t = np.concatenate(([bin_n_t[0]],bin_n_t))
	t_p = bin_n_t/t.size

	a_edges = np.arange(a.min(),a.max()+0.001,0.001)
	bin_n_a,bin_edge_a = np.histogram(a,bins=a_edges)
	bin_edge_a_c = (bin_edge_a[1:]+bin_edge_a[:-1])*0.5
	a_local,a_err_len = get_sigma(bin_edge_a_c,bin_n_a/bin_n_a.sum(),standard=0.96)
	a_start = a_local-a_err_len
	a_stop = a_local+a_err_len
	a_local1,a_err_len1 = get_sigma(bin_edge_a_c,bin_n_a/bin_n_a.sum(),standard=0.68)
	a_err_l = p[1]-a_local1+a_err_len1
	a_err_h = a_local1+a_err_len1-p[1]
	a_edges = np.linspace(a_start,a_stop,50)
	bin_n_a,bin_edge_a = np.histogram(a,bins=a_edges)
	bin_n_a = np.concatenate(([bin_n_a[0]],bin_n_a))
	a_p = bin_n_a/a.size

	eq_edges = np.arange(logeq.min(),logeq.max()+0.01,0.01)
	bin_n_eq,bin_edge_eq = np.histogram(logeq,bins=eq_edges)
	bin_edge_e_c = (bin_edge_eq[1:]+bin_edge_eq[:-1])*0.5
	e_local,e_err_len = get_sigma(bin_edge_e_c,bin_n_eq/bin_n_eq.sum(),standard=0.96)
	e_start = e_local-e_err_len
	e_stop = e_local+e_err_len
	e_local1,e_err_len1 = get_sigma(bin_edge_e_c,bin_n_eq/bin_n_eq.sum(),standard=0.68)
	e_err_l = np.log10(p[2])-e_local1+e_err_len1
	e_err_h = e_local1+e_err_len1-np.log10(p[2])
	eq_edges = np.linspace(e_start,e_stop,50)
	bin_n_eq,bin_edge_eq = np.histogram(logeq,bins=eq_edges)
	bin_n_eq = np.concatenate(([bin_n_eq[0]],bin_n_eq))
	eq_p = bin_n_eq/logeq.size

	c = {}
	c['function'] = function1
	c['p_array'] = [t,a,logeq]
	c['p'] = p
	c['t'] = [p[0],[[t_err_l],[t_err_h]]]
	c['a'] = [p[1],[[a_err_l],[a_err_h]]]
	c['logeq'] = [np.log10(p[2]),[[e_err_l],[e_err_h]]]
	c['t_bound'] = [t_start,t_stop]
	c['a_bound'] = [a_start,a_stop]
	c['e_bound'] = [e_start,e_stop]
	c['hist_t'] = [bin_edge_t,t_p]
	c['hist_a'] = [bin_edge_a,a_p]
	c['hist_eq'] = [bin_edge_eq,eq_p]
	return c




def gbm_lag_analysis(ni_time,ni_energy,bgo_time,bgo_energy,window = None,time_edges = None,dt = 0.1,nai_energy_edges = None,bgo_energy_edges = None,savedir = None):
	if(nai_energy_edges is None):
		nai_energy_edges = np.array([8,10,12,16,21,26,34,43,55,70,90,115,146,186,238,303,386,492,627,800])
	if(bgo_energy_edges is None):
		bgo_energy_edges = np.array([800,1143,1635,2339,3345,4783,6839,9780,13986,20000])
	if(window is None):
		print('window using default method.')
		ni_start,ni_stop = get_pulse_duration(ni_time)
		if(len(ni_start) == 0):
			print('we can not find pulse in this data.')
			return []
		else:
			window = [ni_start[0]-20,ni_stop[-1]+20]
	if(time_edges is None):
		time_edges = [ni_time.min(),ni_time.max()]
	print('window :',window)
	print('dt:',dt)
	print('time edges:',time_edges)
	binarray = np.arange(time_edges[0],time_edges[1]+dt,dt)

	n1 = len(nai_energy_edges)-1
	n2 = len(bgo_energy_edges)-1
	e_index = np.where((ni_energy>=nai_energy_edges[0])&(ni_energy<=nai_energy_edges[1]))
	t_e = ni_time[e_index]
	bin_n,bin_edges = np.histogram(t_e,bins = binarray)
	bin_t = (bin_edges[1:]+bin_edges[:-1])*0.5
	wind_index = np.where((bin_t>=window[0])&(bin_t<=window[1]))
	bin_err_a = np.sqrt(bin_n)
	t_b_a,cs_a,bs = r_baseline(bin_t,bin_n)
	E0 = np.sqrt(nai_energy_edges[1]*nai_energy_edges[0])


	lag = []
	lag_err = []
	egy = []
	for i in range(1,n1):

		e_index = np.where((ni_energy>=nai_energy_edges[i])&(ni_energy<=nai_energy_edges[i+1]))
		t_e= ni_time[e_index]
		bin_n,bin_edges = np.histogram(t_e,bins = binarray)
		bin_t = (bin_edges[1:]+bin_edges[-1])/2.
		bin_err = np.sqrt(bin_n)
		t_b,cs,bs = r_baseline(bin_t,bin_n)
		lag_,err_ = get_lag(cs_a[wind_index],cs[wind_index],dt,x1_err=bin_err_a[wind_index],x2_err=bin_err[wind_index],name= 'n'+str(i),savedir=savedir)

		lag.append(lag_)
		lag_err.append(err_)
		eg = np.sqrt(nai_energy_edges[i+1]*nai_energy_edges[i])
		egy.append(eg)

	cs_a = cs
	bin_err_a = bin_err
	lag_a = lag_
	err_a = err_
	for i in range(n2):

		e_index = np.where((bgo_energy>=bgo_energy_edges[i])&(bgo_energy<=bgo_energy_edges[i+1]))
		t_e = bgo_time[e_index]
		bin_n,bin_edges = np.histogram(t_e,bins = binarray)
		bin_t = (bin_edges[1:]+bin_edges[-1])/2.
		bin_err = np.sqrt(bin_n)
		t_b,cs,bs = r_baseline(bin_t,bin_n)
		lag_, err_ = get_lag(cs_a[wind_index], cs[wind_index],dt, bin_err_a[wind_index] , bin_err[wind_index], name= 'b'+str(i),savedir=savedir)
		lag.append(lag_+ lag_a)
		lag_err.append(err_+err_a)
		eg = np.sqrt(bgo_energy_edges[i+1]*bgo_energy_edges[i])
		egy.append(eg)
	c = {}
	c['lag'] = [np.array(egy),np.array(lag),np.array(lag_err)]
	c['E0'] = E0
	return c










