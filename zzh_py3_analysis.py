import numpy as np
from scipy import optimize


def find_baselin(x,y):
    #we do the frist time
    x0,y0 = find_down(x,y)
    x1,y1 = find_up(x,y)
    #we do zhe second time
    x0,y0 = find_down(x0,y0)
    x1,y1 = find_down(x1,y1)
    #we do zhe second time
    x0,y0 = find_down(x0,y0)
    x1,y1 = find_down(x1,y1)
    # ---------------------
    x_all = np.concatenate((x0, x1))
    y_all = np.concatenate((y0, y1))
    #--------------------
    index_all = np.argsort(x_all)
    x_all = x_all[index_all]
    y_all = y_all[index_all]
    print(x_all)
    print(y_all)
    #print(optimize.curve_fit(f_e,x_all,y_all)[0])
    a,b,c,d,e = optimize.curve_fit(f_4,x_all,y_all)[0]
    #print(a,b,c,d,e)
    bace = [f_4(i,a,b,c,d,e) for i in x]
    bace = np.array(bace)
    y_rate = y - bace
    return x,y,y_rate,bace

def f_4(x,A,B,C,D,E):
    return A*x**4+B*x**3+C*x**2+D*x+E

def find_flux(t,bin,step,t_start,t_stop):
    #t_good = t[np.where(t>=t_start & t<=t_stop)]
    x = np.arange(t_start,t_stop,step)
    y = []
    for i in x:
        k = len(t[np.where((t>=i-bin/2) & (t<=i+bin/2))])/bin
        y.append(k)
    y = np.array(y)
    #print('----')
    #print(len(x))
    #print(len(y))
    return x,y


def find_up(x,y):
    n = len(x)-1
    x0 = []
    y0 = []
    for i in range(1,n):
        if((y[i] > y[i-1]) & (y[i] > y[i+1])):
            x0.append(x[i])
            y0.append(y[i])
    x0 = np.array(x0)
    y0 = np.array(y0)
    return x0,y0

def find_down(x,y):
    n = len(x)-1
    x0 = []
    y0 = []
    for i in range(1,n):
        if((y[i] < y[i-1]) & (y[i] < y[i+1])):
            x0.append(x[i])
            y0.append(y[i])
    x0 = np.array(x0)
    y0 = np.array(y0)
    return x0,y0

