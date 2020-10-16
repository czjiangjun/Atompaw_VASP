#!/usr/bin/env python   
#!coding = utf-8
import os
import time
import PyGnuplot  

# g = PyGnuplot.Gnuplot()   

import PyGnuplot as gp 
import numpy as np 

# X=np.arange(10)
# Y=np.sin(X/(2*np.pi))
# Z=Y**2.0

# gp.s([X,Y,Z])

#gp.c('plot "tmp.dat" u 1:2 w lp')
#gp.c('replot "tmp.dat" u 1:3 w lp')
# gp.c('set encoding utf8')
# gp.c('set xlabel "横坐标"')


# current_path = os.system("pwd")
# print(current_path)
def Write_Atom_in(E_name, Z_Val):
      file_handle=open('in',mode='w+')
      element = [E_name, ' ', str(Z_Val)]
      Atom = ''.join(element)
      file_handle.write(Atom)
      AtomPAW_In = ''' 
PBE-PW  scalarrelativistic loggridv4  2001 5.07
3 3 0 0 0 0
3 1 2
0 0 0
c
c
v
c
v
1
1.85  1.3   1.5    1.5
y
-0.50
n
y
0.22
n
vasprrkj  gramschmidtortho  Besselshape
1 2   MTROULLIER
1.3
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

def Data_Fig(Properties, type, data, name):
    TYPE = type+"-"+data
#    TYPE2 = "ATOM-"+type

    plot_command = "plot "+'"'+type+'_'+data +'" '+'u 1:2 w p pt 3 ps 1.2 title "'+TYPE +'"'
#    plot_command_ATOM = "replot "+'"ATOM_'+data +'" '+'u 1:2 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + data +'"'
    Fig_name = type+"_"+name+".eps"

    gp.c(Fig_title)
    if (Properties == "wavefunction"):
       gp.c('set xrange[0:1.5]')
       gp.c('set yrange[-1.5:1.5]')

    gp.c(plot_command)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
    gp.p(Fig_name)

def Data_compare_Fig(Properties, data, name):
    TYPE1 = "VASP-"+data
    TYPE2 = "ATOM-"+data

    plot_command_VASP = "plot "+'"VASP_'+data +'" '+'u 1:2 w p pt 3 lw 1.2 title "'+TYPE1+'"'
    plot_command_ATOM = "replot "+'"ATOM_'+data +'" '+'u 1:2 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + data +'"'
    Fig_name = name+".eps"

    gp.c(Fig_title)
    if (Properties == "wavefunction"):
       gp.c('set xrange[0:1.5]')
       gp.c('set yrange[-1.5:1.5]')
    gp.c('set yrange[-50:100]')
    gp.c(plot_command_VASP)
    gp.c(plot_command_ATOM)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
    gp.p(Fig_name)

def main():
    Write_Atom_in('Si', 14)

    os.system("~/Softs/atompaw_test/src/atompaw < in > out")
    #exit ()
#    time.sleep(5)

#    Data_Fig("potential", "VASP", "POTAE", "POTAE_data")
#    Data_compare_Fig("density", "CORE", "CORE_data")
#    Data_compare_Fig("density", "PCORE", "PCORE_data")
#    Data_compare_Fig("wavefunction", "WAE", "WAE_data")
#    Data_compare_Fig("wavefunction", "WPS", "WPS_data")
    Data_compare_Fig("potential", "POTAE", "POTAE_data")

if __name__ == '__main__':
    main()
    # print(__name__)

