from astropy.io import fits
import numpy as np
import os
#from zjh_data_analysis import r_baseline
from Separate_background.background_kernel import Baseline_in_time

from glob import glob
import matplotlib.pyplot as plt
import shutil

def make_phaI(bn,ni,topdir,savedir,slice_start,slice_stop,binsize = 1,time_start = -100,time_stop=300):
	'''

	:param bn: 样本代号
	:param ni: 探头
	:param topdir: 数据所在目录
	:param savedir: 产生数据存放目录
	:param slice_start: 光谱切片起始时间
	:param slice_stop: 光谱切片结束时间
	:param binsize: 光变曲线切片大小
	:param time_start: 总数据裁剪，开始时间
	:param time_stop: 总数据裁剪，结束时间

	'''
	datalink = topdir+bn + '/'+'glg_tte_'+ni+'_' + bn + '_v*.fit'
	datalink = glob(datalink)[0]
	if(os.path.exists(savedir) == False):
		os.makedirs(savedir)
	rsp_link = topdir+bn+'/'+'glg_cspec_'+ni+'_'+bn+'_*.rsp'
	rsp_link = glob(rsp_link)[0]
	copy_rspI(rsp_link,savedir+'A_'+bn+'_'+ni+'.rsp')
	hdu = fits.open(datalink)
	trigtime = hdu['Primary'].header['TRIGTIME']
	data_ = hdu['EVENTS'].data

	t = data_.field(0) - trigtime
	ch = data_.field(1)

	ebound = hdu['EBOUNDS'].data

	#ch_n = ebound.field(0)
	emin = ebound.field(1)
	emax = ebound.field(2)

	e_diff = emax-emin
	'''
	首先需要分能道扣除背景，然后设置能谱切片时间
	在这里，设置切片时间是通过设置开始时间点和结束时间点来控制切片。
	
	先获得各个能道的光变曲线，然后通过r_baseline获得背景曲线，最后总计切片中的内容
	
	'''

	edges = np.arange(time_start,time_stop+binsize,binsize)    #生成统计时间片边界数组

	total_rate = np.zeros(128)
	bkg_rate = np.zeros(128)

	total_uncertainty = np.zeros(128)
	bkg_uncertainty = np.zeros(128)

	exposure = len(np.where((edges[:-1]>=slice_start) & (edges[:-1]<=slice_stop))[0][:-1])*binsize

	for i in range(128):
		t_ch = t[ch == i]
		bin_n,bin_edges = np.histogram(t_ch,bins = edges)

		#bin_t = (bin_edges[1:]+bin_edges[:-1])*0.5          #获得单能道光变曲线
		bin_t = bin_edges[:-1]                               #只要前边界
		#bin_rate = bin_n/binsize

		#t_r,cs,bs = r_baseline(bin_t,bin_n)                 #获得单能道背景
		t_r,cs,bs = Baseline_in_time(bin_t,bin_n).get_value()
		slice_index = np.where((bin_t>=slice_start) & (bin_t<=slice_stop))[0]
		slice_index = slice_index[:-1]                      #排除掉最后一个可能不完整的bin
		#print('slice index:\n',slice_index)
		#print('bs:\n',bs[slice_index])
		total_rate[i] = (bin_n[slice_index]).sum()#.mean()
		bkg_rate[i] = (bs[slice_index]).sum()#.mean()
		if(total_rate[i] <= bkg_rate[i]):
			bkg_rate[i] = total_rate[i] #限制背景高度

		#exposure = len(slice_index)*binsize
		bkg_uncertainty[i] = np.sqrt(bkg_rate[i]/exposure)
		total_uncertainty[i] = np.sqrt(total_rate[i]/exposure)

	write_phaI(total_rate,bn,ni,slice_start,slice_stop,savedir+'A_'+bn+'_'+ni+'.pha', exposure)
	write_phaI(bkg_rate,bn,ni,slice_start,slice_stop,savedir+'A_'+bn+'_'+ni+'.bkg',exposure)

	x = np.sqrt(emin*emax)

	plt.figure(figsize=(10,10))
	plt.subplot(1,1,1)
	plt.errorbar(x,bkg_rate/e_diff,yerr = bkg_uncertainty/e_diff,color = 'blue')
	plt.errorbar(x,total_rate/e_diff,yerr = total_uncertainty/e_diff,color = 'r')
	plt.xlabel('energy KeV')
	plt.ylabel('counts /N')
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(savedir+'Z_slic_'+bn+'_'+ni+'.png')
	plt.close()


def write_phaI(spectrum_counts,bnname,detector,t1,t2,outfile,exposure):
	'''

	:param spectrum_counts: 光谱计数
	:param bnname: 样本代号
	:param detector: 探头
	:param t1: 光谱切片起始时间
	:param t2: 光谱切片结束时间
	:param outfile: 输出文件保存路径
	:param exposure: 光谱曝光时间
	'''
	header0=fits.Header()#头文件基本信息设置
	header0.append(('creator', 'Zou', 'The name who created this PHA file'))
	header0.append(('telescop', 'Fermi', 'Name of mission/satellite'))
	header0.append(('bnname', bnname, 'Burst Name'))
	header0.append(('t1', t1, 'Start time of the PHA slice'))
	header0.append(('t2', t2, 'End time of the PHA slice'))

	hdu0=fits.PrimaryHDU(header=header0) #创建头

	a1 = np.arange(128)
	col1 = fits.Column(name='CHANNEL', format='1I', array=a1)                              #创建列
	col2 = fits.Column(name='COUNTS', format='1D', unit='COUNTS', array=spectrum_counts)     #创建列
	hdu1 = fits.BinTableHDU.from_columns([col1, col2])#创建一个bin列表
	header=hdu1.header
	header.append(('extname', 'SPECTRUM', 'Name of this binary table extension'))
	header.append(('telescop', 'GLAST', 'Name of mission/satellite'))
	header.append(('instrume', 'GBM', 'Specific instrument used for observation'))
	header.append(('filter', 'None', 'The instrument filter in use (if any)'))
	header.append(('exposure', exposure, 'Integration time in seconds'))#这里设置曝光时间
	header.append(('areascal', 1., 'Area scaling factor'))
	header.append(('backscal', 1., 'Background scaling factor'))
	if outfile[-3:]=='pha':
		header.append(('backfile', 'A_'+bnname+'_'+detector+'.bkg', 'Name of corresponding background file (if any)'))#这里有背景文件的名字
		header.append(('respfile', 'A_'+bnname+'_'+detector+'.rsp', 'Name of corresponding RMF file (if any)'))#这里有响应文件的名字
	else:
		header.append(('backfile', 'none', 'Name of corresponding background file (if any)'))
		header.append(('respfile', 'none', 'Name of corresponding RMF file (if any)'))
	header.append(('corrfile', 'none', 'Name of corresponding correction file (if any)'))
	header.append(('corrscal', 1., 'Correction scaling file'))
	header.append(('ancrfile', 'none', 'Name of corresponding ARF file (if any)'))
	header.append(('hduclass', 'OGIP', 'Format conforms to OGIP standard'))
	header.append(('hduclas1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'))
	header.append(('hduclas2', 'TOTAL', 'Indicates gross data (source + background)'))
	header.append(('hduclas3', 'COUNT', 'Indicates data stored as counts'))
	header.append(('hduvers', '1.2.1', 'Version of HDUCLAS1 format in use'))
	header.append(('poisserr', True, 'Use Poisson Errors (T) or use STAT_ERR (F)'))
	header.append(('chantype', 'PHA', 'No corrections have been applied'))
	header.append(('detchans', 128, 'Total number of channels in each rate'))
	header.append(('hduclas4', 'TYPEI', 'PHA Type I (single) or II (mulitple spectra)'))

	header.comments['TTYPE1']='Label for column 1'
	header.comments['TFORM1']='2-byte INTERGER'
	header.comments['TTYPE2']='Label for column 2'
	header.comments['TFORM2']='8-byte DOUBLE'
	header.comments['TUNIT2']='Unit for colum 2'

	hdul = fits.HDUList([hdu0, hdu1])#保存
	if(os.path.exists(outfile)):
		os.remove(outfile)#删除旧版本
	hdul.writeto(outfile)
	hdul.close()

def copy_rspI(datalink,outfile):

	if(os.path.exists(outfile)): #检查文件是否存在，如果存在删除旧版
		os.remove(outfile)
	shutil.copyfile(datalink,outfile)







