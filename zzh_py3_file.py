import re
#import numpy as np
import os


def transpose(data):
    da = [[row[col] for row in data] for col in range(len(data[0]))]
    return da

def printdatatofile(file_name , data):
    data_transpose = transpose(data)
    f = open(file_name,'w+')
    number = len(data[0])
    for j in range(number):
        data_h = data_transpose[j]
        n = len(data_h)
        strne = '\t'+str(data_h[0])
        for k in range(1,n):
            strne = strne + '\t' + str(data_h[k])
        strne = strne + '\n'
        f.write(strne)
    f.close()

def readcol(file_name):
    c = []
    f = open(file_name,'r')
    for line in f:
        a = line
        b = [i for i in re.split('[\t,\n,\s]',a) if i != '']
        #print(b)
        c = c + [b]
    f.close()
    c = [i for i in c if i != []]
    data_all = transpose(c)
    nl = len(data_all)
    for k in range(nl):
        for l in range(len(data_all[k])):
            try:
                data_all[k][l] = int(data_all[k][l])
            except ValueError:
                try:
                    data_all[k][l] = float(data_all[k][l])
                except ValueError:
                    data_all[k][l] = data_all[k][l]
    return data_all

def getfilelist(dir1):
    if (os.path.exists(dir1)):
        dirnames = os.listdir(dir1)
        dirlist = []
        dir_number = 0
        filelist = []
        file_number = 0
        for sample in dirnames:
            if (os.path.isfile(dir1 + sample)):
                filelist.append(sample)
                file_number = file_number + 1
            elif(os.path.isdir(dir1 + sample)):
                dirlist.append(sample)
                dir_number = dir_number + 1
        c = {dir1: filelist}
        if (dir_number != 0):
            for sname in dirlist:
                dir2 = dir1+sname+'/'
                k = getfilelist(dir2)
                c.update(k)
        return c
    else:
        print('do not find the dir named [' + dir1 + ']\n')
        return False

def findfile(dir1,feature):
    if (os.path.exists(dir1)):
        dirnames = os.listdir(dir1)
        filelist = []
        fil_number = 0
        fil_result_number = 0
        featurelist = [i for i in re.split('[*]',feature) if i != '']
        #print(featurelist)
        for_number = len(featurelist)
        fileresult = [[] for i in range(for_number)]
        for eve in range(for_number):
            if(eve == 0):

                fe_number = len(featurelist[eve])
                for sample in dirnames:
                    if (os.path.isfile(dir1 + sample)):
                        filelist.append(sample)
                        fil_number = fil_number + 1
                if (fil_number != 0):
                    for i in filelist:
                        i_number = len(i)
                        n = i_number - fe_number + 1
                        for j in range(n):
                            if (i[j:j + fe_number] == featurelist[eve]):
                                fileresult[eve].append(i)
                                fil_result_number = fil_result_number + 1
                                break
                    #print('1----------',fileresult[eve])#------------------------
                    if (fil_result_number == 0):
                        print('we do not find any file that has the feature with [' + feature + ']!\n')
                        return (0)
                    else:
                        fil_result_number = 0
                else:
                    print('there is no file in this dir ! \n')
                    return 0
            else:
                fe_number = len(featurelist[eve])
                for i in fileresult[eve-1]:
                    i_number = len(i)
                    n = i_number - fe_number + 1
                    for j in range(n):
                        if (i[j:j + fe_number] == featurelist[eve]):
                            fileresult[eve].append(i)
                            fil_result_number = fil_result_number + 1
                            break
                #print('2----------', fileresult[eve])  # ------------------------
                if (fil_result_number == 0):
                    print('we do not find any file that has the feature with [' + feature + ']!\n')
                    return (0)
                else:
                    fil_result_number = 0
        return fileresult[for_number-1]
    else:
        print('do not find the dir named ['+ dir1 + ']!\n')
        return False

#dir = '/home/laojin/my_work/'
#name = 'my_sample.txt'
#resultdata = readcol(dir+name)

#a = findfile(dir,'tte*_n0')
#print(resultdata)




