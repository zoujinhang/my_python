import os
from ftplib import FTP_TLS as FTP
from multiprocessing import Pool


def download_all_in_one_path(targetdir,resultdir,check = True,num = 50):
	if(os.path.exists(resultdir) == False):
		os.makedirs(resultdir)
	ftp = FTP('129.164.179.23')
	ftp.login()
	ftp.prot_p()
	ftp.cwd(targetdir)
	files = ftp.nlst()
	target = 'https://heasarc.gsfc.nasa.gov/FTP' + targetdir
	c = None
	if(check):
		c = []
	data1 = []
	ftp.voidcmd('TYPE I')
	print('正在获取校验信息........')
	for i in files:
		#print(i)
		data = os.path.join(target,i)
		print(data)
		data1.append(data)
		if(check):
			c.append(ftp.size(i))
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
			for index,i in enumerate(targlist):
				nn = True
				for j in eee:
					if (os.path.split(i)[1] == j):
						if ((check == None) | (targnumber != len(check))):
							nn = False
							break
						else:
							myfilesize = os.path.getsize(j)
							if (myfilesize >= check[index]):
								nn = False
								break
							else:
								os.system('rm '+j)
								break
				if (nn):
					kkk = check[index]
					en.append([i,kkk,resultdir])
		targlist = en
	else:
		en = []
		for index,i in enumerate(targlist):
			kkk = check[index]
			en.append([i,kkk,resultdir])
		targlist = en
	print(targlist)
	pool = Pool(threadnum)
	pool.map(download,targlist)
	pool.close()
	pool.join()
	print('目标下需要载数：',targnumber)
	print('重复目标数：',targnumber-len(targlist))
	print('实际下需要载数：',len(targlist))


def download(target1):
	target = target1[0]
	check = target1[1]
	local = target1[2]
	filename = os.path.split(target)[1]
	t_link = 'wget --quiet --show-progress --read-timeout=5 --tries=0 -P '+local + ' '+target
	#print(t_link)
	os.system(t_link)
	if(os.path.exists(filename) == False):
		print('\n' + filename + ' 下载失败，即将重新下载！')
		download(target1)
	else:
		myfilesize = os.path.getsize(filename)
		if(myfilesize < check):
			print('\n'+filename + '未通过校验，即将重新下载！！')
			os.system('rm '+filename)
			download(target1)

#targetdir = '/fermi/data/gbm/daily/2017/02/08/current/'
#resultdir = '/home/laojin/gbm_daily_database/2017/02/08/'
#download_all_in_one_path(targetdir,resultdir,num = 10)





































