from xspec import *
import os
import matplotlib.pyplot as plt
import numpy as np
def xspec_fit_kernel(filelist,datadir,savedir):
	os.chdir(datadir)
	alldatastr = ' '.join(filelist)
	AllData.clear()
	AllModels.clear()
	AllData(alldatastr)
	AllData.show()
	AllData.ignore('1:**-200.0,40000.0-** 2-3:**-8.0,800.0-**')
	Model('grbm')
	Fit.statMethod='pgstat'
	Fit.nIterations=1000
	Fit.query = "yes"
	Fit.perform()
	Fit.error('3.0 3')
	Fit.perform()
	AllModels.calcFlux("8. 40000.0 err") #参数需要一个能量的范围，之前在拟合光谱时我们设置了一个范围，我们暂时用它。
	par3=AllModels(1)(3)#第一个模型的第三个参数
	value = par3.values[0]
	value_arr1,value_arr2,ffff = par3.error
	Plot('eeufspec')
	flux_list = []
	for i in range(len(filelist)):
		print(i)
		flux = AllData(i+1).flux
		flux_list.append(flux)
		energies=Plot.x(i+1)
		rates=Plot.y(i+1)
		folded=Plot.model(i+1)
		xErrs=Plot.xErr(i+1)
		yErrs=Plot.yErr(i+1)
		plt.errorbar(energies,rates,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
		plt.plot(energies,folded,color='black',zorder=2)
	plt.axvline(x = value,color = 'r')
	plt.axvline(x = value_arr1,color = 'g')
	plt.axvline(x = value_arr2,color = 'g')
	plt.xlabel('Energy KeV')
	plt.ylabel(r'${KeV^{2} (Photons cm^{-2}s^{-1}keV^{-1})}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(savedir + 'foldedspec.png')
	#plt.close()
	return value,value_arr1,value_arr2,np.array(flux_list).T
















