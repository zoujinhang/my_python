import os
import queue
import threading
from ftplib import FTP
#from multiprocessing import Pool

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

    print('目标下需要载数：',targnumber)
    print('重复目标数：',targnumber-len(targlist))
    print('实际下需要载数：',len(targlist))

    qu = queue.Queue()
    for i in targlist:
        print(i)
        qu.put(i)
    num = threading.Semaphore(threadnum)
    n = len(targlist)
    m = 0
    thread = []
    for i in range(n):
        if((check == None) | (targnumber != len(check))):
            t = Download(qu,num)
        else:
            t = Download(qu,num,check)
        t.start()
        thread.append(t)
    for i in thread:
        i.join()
        m = m+1
    #print('hehehehehehehehehhehehehehehehehehehehe!!')
    qu.join()
    print('\n下载完成数：',m)


class Download(threading.Thread):
    def __init__(self,list,num,check = None):
        super().__init__()
        self.list = list
        self.num = num
        self.check = check
    def run(self):
        with self.num:
            targ = self.list.get()
            filename = os.path.split(targ)[1]
            print('\n线程'+self.name)
            print('\n目标下载以开始：'+ filename)
            download(targ,self.check[filename])
            print('\n'+filename + '下载完成！！')
            self.list.task_done()


def download(target,check = None):
    filename = os.path.split(target)[1]
    t_link = 'wget --quiet --show-progress --read-timeout=5 --tries=0 '+target
    os.system(t_link)
    if(os.path.exists(filename) == False):
        print('\n' + filename + ' 下载失败，即将重新下载！')
        download(target,check)
    else:
        if(check != None):
            myfilesize = os.path.getsize(filename)
            if(myfilesize < check):
                print('\n'+filename + '未通过校验，即将重新下载！！')
                print('\n实际文件大小：',check)
                print('\n下载文件大小：',myfilesize)
                os.system('rm '+filename)
                download(target,check)


#targetdir = '/fermi/data/gbm/bursts/2018/bn180426549/current/'
#resultdir = '/home/laojin/my_work/data/bn180426549/'
#download_all_in_one_path(targetdir,resultdir)


