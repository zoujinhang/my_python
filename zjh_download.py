import os
from ftplib import FTP
from multiprocessing import Pool


def download_all_in_one_path(targetdir,resultdir,check = True,num = 50):
	if(os.path.exists(resultdir) == False):
		os.makedirs(resultdir)
	ftp = FTP('129.164.179.23')
	ftp.login()
	ftp.cwd(targetdir)
	files = ftp.nlst()
	target = 'https://heasarc.gsfc.nasa.gov/FTP' + targetdir
	c = None
	if(check):
		c = {}
	data1 = []
	ftp.voidcmd('TYPE I')
	print('正在获取校验信息........')
	for i in files:
		print(i)
		data = os.path.join(target,i)
		data1.append(data)
		if(check):
			c[i] = ftp.size(i)
	ftp.quit()
	if(check == False):
		print('忽略数据大小校验。')
	print('正在校验...............')
	down(data1,resultdir,check=c,threadnum = num)
	print('\n任务下载完成！！！')


def down(targlist,resultdir,check = None,threadnum = 50):
	if(os.path.exists(resultdir) == False):
		os.makedirs(resultdir)
	os.chdir(resultdir)
	targnumber = len(targlist)
	rea = os.listdir(resultdir)
	nu = len(rea)
	if (nu != 0):
		eee = []
		for i in rea:
			if (os.path.isfile(i)):
				eee.append(i)
		if (len(eee) != 0):
			en = []
			for i in targlist:
				nn = True
				for j in eee:
					if (os.path.split(i)[1] == j):
						if ((check == None) | (targnumber != len(check))):
							nn = False
							break
						else:
							myfilesize = os.path.getsize(j)
							if (myfilesize >= check[os.path.split(i)[1]]):
								nn = False
								break
							else:
								os.system('rm '+j)
								break
				if (nn):
					en.append(i)
		targlist = en
	pool = Pool(threadnum)
	pool.map(download,targlist)
	print('目标下需要载数：',targnumber)
	print('重复目标数：',targnumber-len(targlist))
	print('实际下需要载数：',len(targlist))


def download(target):
	t_link = 'wget --quiet --show-progress --read-timeout=5 --tries=0 '+target
	os.system(t_link)

#targetdir = '/fermi/data/gbm/bursts/2018/bn180906759/current/'
#resultdir = '/home/laojin/my_work/data/bn180906759/'
#download_all_in_one_path(targetdir,resultdir)





































