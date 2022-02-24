#!/usr/bin/env python   
#!coding = utf-8
import os
import time
import random

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

class Fig_linestyle:
      def __init__(self, style, index=None, line=None):
          self.index = list(style.keys())
          self.line = list(style.values())

class Fig_Draw:
      def __init__(self, name, file, xlabel= None, ylabel= None, cure_type = None, marker= None):
          self.name = name
          self.file = file
          self.xlabel = xlabel
          self.ylabel = ylabel
          self.cure_type = cure_type
          self.marker = marker


      def Data(self, file, index):
          Data_XX = []
          Data_YY = []
          Data_File = open(file, 'r')
          num = 0
          for line in Data_File:
              Data = line.replace('\n', '').strip('').split(' ')
              Data = [i for i in Data if i !='']
#               print(Data[1])
              Data_XX.append(float(Data[0]))
              Data_YY.append(float(Data[index]))
              num = num + 1
          return Data_XX, Data_YY, num

      def Fig(self, Data_XX, Data_YY):
          xmin = round(Data_XX[0])-0.05
          xmax = round(Data_XX[-1])
   
          ymin = round(min(Data_YY))
          ymax = round(max(Data_YY))
          ymin = ymin-(ymax-ymin)/5.0
          ymax = ymax+(ymax-ymin)/5.0

          plt.axis([xmin,xmax,ymin,ymax])
          plt.xlabel(self.xlabel)
          plt.ylabel(self.ylabel)
          plt.plot(Data_XX, Data_YY, label = self.name, ls = self.cure_type, marker = self.marker)
#           plt.legend(loc='best') 
          plt.title(self.name)

#          plt.savefig(local_DOS_fig)
          return

#     
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

def Fig_wave_Draw(Data_XX, Data_YY, title, label_ind, xlabel, ylabel, type, maker):
    xmin = round(Data_XX[0])-0.05
    xmax = round(Data_XX[-1])
    
    ymin = round(min(Data_YY))-0.1
    ymax = round(max(Data_YY))+0.5

    plt.axis([xmin,xmax,ymin,ymax])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(Data_XX, Data_YY, label = label_ind, ls = type, marker = maker)

    plt.legend(loc='best') 
    plt.title(title)

#    plt.savefig(local_DOS_fig)
#    plt.show()
    return

def Fig_potential_Draw(Data_XX, Data_YY, title, label_ind, xlabel, ylabel, type, maker):
    xmin = round(Data_XX[0])-0.05
    xmax = round(Data_XX[-1])
    
    ymin = round(min(Data_YY))
    ymax = round(max(Data_YY))
    ymin = ymin-(ymax-ymin)/10.0
    ymax = ymax+(ymax-ymin)/10.0

    plt.axis([xmin,xmax,ymin,ymax])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(Data_XX, Data_YY, label = label_ind, ls = type, marker = maker)

    plt.legend(loc='best') 
    plt.title(title)

#    plt.savefig(local_DOS_fig)
#    plt.show()
    return

def Get_Fig_Wave_Data(File, types, type, maker):
    Data_XX = []
    Data_YY = []
    tot_lines = os.popen("wc -l " + File).readlines()[0].replace('\n','').strip('').split(' ')
#    print(tot_lines)
    num_unit = int(int(tot_lines[0])/int(types))

    Data_File = open(File, 'r')
    num = 0
#    times = -1
    for line in Data_File:
        num = num + 1
        Data = line.replace('\n', '').strip('').split(' ')
        Data = [i for i in Data if i !='']
#        print(Data)
        mod = num % num_unit
#        print(mod)
        if (mod !=0):
            Data_XX.append(float(Data[0])) 
            Data_YY.append(float(Data[1]))
        else:
#            print(Data_XX, Data_YY)
#            exit(0)
            xlabel = 'Radial /A'
            ylabel = 'Wavefunction'
            Fig_wave_Draw(Data_XX, Data_YY, 'WAE', 'wavefunction', xlabel, ylabel, type, maker)
            Data_XX = []
            Data_YY = []
    return 


def main():
    Write_Atom_in('Si', 14)

    os.system("~/Softs/atompaw_test/src/atompaw < in > out")
    #exit ()

    style = {1:'-',2:'-.',3:'--',4:'-'}
    style_index = list(style.keys())
    style_line = list(style.values())

    type = style_line[random.choice(style_index)-1]
#    print(type)
    Get_Fig_Wave_Data('ATOM_WAE', 4, type, '*')
#    exit(0)
#    type = style_line[random.choice(style_index)-1]
    Get_Fig_Wave_Data('VASP_WAE', 4, '', '+')
#    plt.savefig(local_DOS_fig)
    plt.show()

    FILE_NAME = 'VASP_CORE'
    xlabel = 'Radial /A'
    ylabel = 'Density'
    for i in range(1):
        XX, YY, lines = Get_Fig_Potential_Data(FILE_NAME, i+1)
#        type = style_line[random.choice(style_index)-1]
        Fig_potential_Draw(XX, YY, FILE_NAME, 'Core_density', xlabel, ylabel, '', '+')
    FILE_NAME = 'ATOM_CORE'
    xlabel = 'Radial /A'
    ylabel = 'Density'
    for i in range(1):
        XX, YY, lines = Get_Fig_Potential_Data(FILE_NAME, i+1)
        type = style_line[random.choice(style_index)-1]
        Fig_potential_Draw(XX, YY, FILE_NAME, 'Core_density', xlabel, ylabel, type, '')
#    plt.savefig(local_DOS_fig)
    plt.show()

    FILE_NAME = 'VASP_POTAE'
    xlabel = 'Radial /A'
    ylabel = 'Potential'
    for i in range(2):
        XX, YY, lines = Get_Fig_Potential_Data(FILE_NAME, i+1)
        type = style_line[random.choice(style_index)-1]
        Fig_potential_Draw(XX, YY, FILE_NAME, 'Potential', xlabel, ylabel, type, '*')
#    plt.savefig(local_DOS_fig)
    plt.show()

    print('class_TEST')
    FILE_NAME = 'VASP_POTPS'
    xlabel = 'Radial /A'
    ylabel = 'Potential /Ry'

    linestyle = Fig_linestyle(style)
    type = random.choice(linestyle.line)
#    print(type)
    Potential = Fig_Draw('Potential', 'VASP_POTCAR', xlabel, ylabel, cure_type = type, marker = '*')
#    Potential = Fig_Draw('Potential', 'VASP_POTCAR', xlabel, ylabel, cure_type = type)
    for i in range(2):
        XX, YY, lines = Potential.Data(FILE_NAME, i+1)
        Potential.Fig(XX, YY)
#        print(XX)
    plt.savefig('POTPS_COMPARE_FIG') #    plt.savefig('POTPS_COMPARE_FIG.eps/pdf')
    plt.show()

if __name__ == '__main__':
    main()
    # print(__name__)

