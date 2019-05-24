import zzh_py3_file as zhf
import numpy as np
import os
from decimal import *

def r_baseline(x,y,result = '',lamb = 1,hwi = 30,it = 20,int = 200):
    if (result != ''):
        if(os.path.exists(result) == False):
            os.mkdir(result)
    else:
        result = os.getcwd()
    dir_d = os.path.join(result, 'in_xy.txt')
    data_d = [x, y]
    zhf.printdatatofile(dir_d, data_d)
    x = np.array(x,np.double)
    y = np.array(y,np.double)
    os.system('Rscript  /home/laojin/my_work/my_r/zzh_baseline_r.r ' + dir_d + ' ' + result + ' '+str(lamb) + ' '+ str(hwi) + ' '+ str(it) + ' '+ str(int) + '>>/dev/null')
    if(os.path.exists(os.path.join(result,'back.txt'))):
        back = np.array(zhf.readcol(os.path.join(result,'back.txt'))[0])
        net = y - back
        x = decimal_for_array(x,20)
        y = decimal_for_array(y,20)
        net = decimal_for_array(net,20)
        back = decimal_for_array(back,20)
        zhf.printdatatofile(os.path.join(result, 'xout.txt'), [x])
        zhf.printdatatofile(os.path.join(result, 'net.txt'), [net])
        zhf.printdatatofile(os.path.join(result, 'back.txt'), [back])
        zhf.printdatatofile(os.path.join(result,'summary_t_obs_net_back.txt'),[x,y,net,back])
        return x,y,net,back
    else:
        print('failed!!')
        print('I can not do the following!!')
        return False

def easy_baseline(x,y,result = ''):
    if (result != ''):
        if(os.path.exists(result) == False):
            os.mkdir(result)
    else:
        result = os.getcwd()
    dir_d = os.path.join(result, 'in_xy.txt')
    data_d = [x, y]
    zhf.printdatatofile(dir_d, data_d)
    x = np.array(x,np.double)
    y = np.array(y,np.double)
    os.system('Rscript  /home/laojin/software/idl_my_lib/zbbidl/baseline_r_2015.r ' + dir_d + ' ' + result + '>>/dev/null')
    if(os.path.exists(os.path.join(result,'back.txt'))):
        back = np.array(zhf.readcol(os.path.join(result,'back.txt'))[0])
        net = y - back
        x = decimal_for_array(x,20)
        y = decimal_for_array(y,20)
        net = decimal_for_array(net,20)
        back = decimal_for_array(back,20)
        zhf.printdatatofile(os.path.join(result, 'xout.txt'), [x])
        zhf.printdatatofile(os.path.join(result, 'net.txt'), [net])
        zhf.printdatatofile(os.path.join(result, 'back.txt'), [back])
        zhf.printdatatofile(os.path.join(result,'summary_t_obs_net_back.txt'),[x,y,net,back])
        return x,y,net,back
    else:
        print('failed!!')
        print('I can not do the following!!')
        return False


def decimal_for_array(array,bit):
    r = []
    for i in array:
        n =Decimal(i)
        n = round(n,bit)
        cc = '{:> %(n1)d}' % {'n1':bit+5}
        n = cc.format(n)
        r.append(n)
    r = np.array(r)
    return r



#result = '/home/laojin/shiyan/python/'
#data = result + 'bn130427324.txt'
#resultd = result + 'bn130427324'
#time,flux = zhf.readcol(data)
#time1,flux,net,back = r_baseline(time,flux,resultd,1,30,20,200)



