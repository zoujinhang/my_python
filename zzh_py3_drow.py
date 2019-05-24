import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def easy_scatter(dir,x,y,xlabel,ylabel):
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.scatter(x,y,c = 'k',marker = ',',s = 1)
    fig.savefig(dir)

def easy_plot(dir,x,y,xlabel,ylabel,log = 'False'):
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    n = len(x)
    #ax.set_xticks(np.linspace(0,0.01,4))
    ax.loglog(x,y)
    A,B,C,D,E = optimize.curve_fit(f_4,x,y)[0]
    x0 = np.linspace(x[0],x[n-1],10000)
    y0 = A*x0**4 + B*x0**3 + C*x0**2 + D*x0 + E
    ax.plot(x0,y0,'green')
    ax.scatter(x,y,c = 'y',marker= 'o',s = 20)
    fig.savefig(dir)

def f_4(x,A,B,C,D,E):
    return A*x**4+B*x**3+C*x**2 +D*x + E

def f_3(x,A,B,C,D):
    return A*x**3 + B*x**2 + C*x + D

def f_2(x, A, B, C):
    return A*x*x + B*x + C

def f_1(x, A, B):
    return A*x + B

def f_gauss(x, A, B, C, sigma):
    return A*np.exp(-(x-B)**2/(2*sigma**2)) + C