#!/usr/bin/env python   
#!coding = utf-8
import os
import time

# g = PyGnuplot.Gnuplot()   

import PyGnuplot as gp 
import numpy as np 

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D

# X=np.arange(10)
# Y=np.sin(X/(2*np.pi))
# Z=Y**2.0

# gp.s([X,Y,Z])

#gp.c('plot "tmp.dat" u 1:2 w lp')
#gp.c('replot "tmp.dat" u 1:3 w lp')
# gp.c('set encoding utf8')
# gp.c('set xlabel "x"')


# current_path = os.system("pwd")
# print(current_path)
def Write_Atom_in(E_name, Z_Val):
      file_handle=open('in',mode='w+')
      element = [E_name, ' ', str(Z_Val)]
      Atom = ''.join(element)
      file_handle.write(Atom)
      AtomPAW_In = ''' 
PBE-PW  scalarrelativistic loggridv4  2001  5 1.55 loggderivrange -25. 25. 2500
3 3 0 0 0 0
3 1 2
0 0 0
c
c
v
c
v
1
1.65  1.51   1.64    1.51  # 1.83  1.41   1.07    1.51
y
-0.60    # 能量往大调压制波函数 
n
y
0.42    # 能量往大调压制波函数 
n
vasprrkj  gramschmidtortho  Besselshape
1 3   MTROULLIER
1.5
1.5
1.6
1.7
1.6
XMLOUT
default
PWPAWOUT
ABINITOUT
default
PWSCFOUT
UPFDX  0.0125d0   UPFXMIN   -7.d0     UPFZMESH 11.d0
END\n  '''
      file_handle.writelines(AtomPAW_In)
      file_handle.flush()
      file_handle.close()

def Get_Fig_Wave_Data(File, num_unit):
    Data_XX = []
    Data_YY = []
    Data_File = open(File, 'r')
    num = 0
    times = -1
    for line in Data_File:
        Data = line.replace('\n', '').strip('').split(' ')
#        print(Data)
        Data = [i for i in Data if i !='']
        mod = num % num_unit
#        print(num)
        if mod != (num_unit -1):
            Data_XX.append(float(Data[0])) 
            Data_YY.append(float(Data[1]))
#        if (num == num_unit):
#            break
#        print(Data_XX, Data_YY)
        num = num + 1
    return Data_XX, Data_YY

def Get_Fig_Potential_Data(File, index):
    Data_XX = []
    Data_YY = []
    Data_File = open(File, 'r')
    num = 0
    for line in Data_File:
        Data = line.replace('\n', '').strip('').split(' ')
        Data = [i for i in Data if i !='']
#        print(Data[1])
        Data_XX.append(float(Data[0]))
        Data_YY.append(float(Data[index]))
        num = num + 1
    return Data_XX, Data_YY, num

def Fig_Draw(Data_XX, Data_YY, title, label_ind, xlabel, ylabel):
    xmin = round(Data_XX[0])-0.2
    xmax = round(Data_XX[-1])
    
    if (Data_YY[0] < Data_YY[-1]):
       ymin = round(Data_YY[0])-1
       ymax = round(Data_YY[-1])+1
    else:
       ymin = round(Data_YY[-1])-1
       ymax = round(Data_YY[0])+1

    plt.axis([xmin,xmax,ymin,ymax])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(Data_XX, Data_YY, label = label_ind)

    plt.legend(loc='best') 
    plt.title(title)

#    plt.savefig(local_DOS_fig)
    plt.show()
    return


def main():
    Write_Atom_in('Si', 14)

    os.system("~/Softs/atompaw_test/src/atompaw < in > out")
    #exit ()

    FILE_NAME = 'VASP_POTAE'
    xlabel = 'Radial /A'
    ylabel = 'Potential'
    for i in range(3):
       XX, YY, lines = Get_Fig_Potential_Data(FILE_NAME, i)
       Fig_Draw(XX, YY, FILE_NAME, 'Potential', xlabel, ylabel)

    XX, YY = Get_Fig_Wave_Data('ATOM_WAE', 2002)
    xlabel = 'Radial /A'
    ylabel = 'Wavefunction'
    Fig_Draw(XX, YY, 'WAE', 'wavefunction', xlabel, ylabel)

#    plt.savefig(local_DOS_fig)
    plt.show()


if __name__ == '__main__':
    main()
    # print(__name__)

