#!/usr/bin/env python   
#!coding = utf-8
import os
import time

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
# gp.c('set xlabel "x"')


# current_path = os.system("pwd")
# print(current_path)
def Write_Atom_in(E_name, Z_Val):
      file_handle=open('in',mode='w+')
      element = [E_name, ' ', str(Z_Val)]
      Atom = ''.join(element)
      file_handle.write(Atom)
      AtomPAW_In = ''' 
LDA-PW scalarrelativistic loggridv4 2001 5 1.55 loggderivrange -25. 25. 2500
3 3 0 0 0 0
3 1 2
0 0 0
c
c
v
c
v
1
1.65  1.21   1.07    1.21
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

def Data_Fig(Properties, type, data, name):
    TYPE = type+"-"+data
#    TYPE2 = "ATOM-"+type

    plot_command1 = "plot "+'"'+type+'_'+data +'" '+'u 1:2 w p pt 3 ps 1.2 title "'+TYPE +'"'
    plot_command2 = "replot "+'"'+type+'_'+data +'" '+'u 1:3 w l lt 5 lw 1.2 title "'+TYPE +'"'
#    plot_command_ATOM = "replot "+'"ATOM_'+data +'" '+'u 1:2 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + data +'"'
    Fig_name = type+"_"+name+".eps"

    gp.c(Fig_title)
    if (Properties == "wavefunction"):
       gp.c('set xrange[0:1.5]')
       gp.c('set yrange[-1.5:1.5]')

    gp.c(plot_command1)
    gp.c(plot_command2)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
    gp.p(Fig_name)

def Data_compare_Fig(Properties, data, name):
    TYPE1 = "VASP-"+data
    TYPE2 = "ATOM-"+data

    plot_command_VASP = "plot "+'"VASP_'+data +'" '+'u 1:2 w p pt 1 lw 1.0 title "'+TYPE1+'"'
    plot_command_ATOM = "replot "+'"ATOM_'+data +'" '+'u 1:2 w l lt 3 lw 1.5 title "'+TYPE2+'"'
#    plot_command_ATOM = "replot "+'"VASP_'+data +'" '+'u 1:3 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + data +'"'
    Fig_name = name+".eps"

    gp.c(Fig_title)
    gp.c('set xrange[0:3.0]')
    if (Properties == "waefunction"):
       gp.c('set yrange[-1:1.5]')
    elif (Properties == "wpsfunction"):
       gp.c('set yrange[-0.5:1.5]')
    elif (Properties == "density" and data == "PCORE"):
       gp.c('set yrange[0:1]')
    elif (Properties == "density" and data == "CORE"):
       gp.c('set yrange[0:9]')
    elif (Properties == "potential" and data == "POTAE"):
       gp.c('set yrange[-50: 5000]')
    elif (Properties == "potential" and data == "POTPS"):
       gp.c('set yrange[-50:20]')
    elif (Properties == "potential" and data == "POTPSC"):
       gp.c('set yrange[-500:100]')
    else:
#       gp.c('set yrange[-50:250]')
       gp.c('set yrange[0:4]')
    gp.c(plot_command_VASP)
    gp.c(plot_command_ATOM)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
    gp.p(Fig_name)

def Data_compare_Fig2(Properties, data1, data2, name):
    TYPE1 = "VASP-"+data1
    TYPE2 = "VASP-"+data2

    plot_command1 = "plot "+'"VASP_'+data1 +'" '+'u 1:2 w p pt 3 lw 1.2 title "'+TYPE1+'"'
    plot_command2 = "replot "+'"VASP_'+data2 +'" '+'u 1:2 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    plot_command_ATOM = "replot "+'"VASP_'+data +'" '+'u 1:3 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + name +'"'
    Fig_name = name+".eps"

    gp.c(Fig_title)
    gp.c('set xrange[0:3.0]')
    if (Properties == "waefunction"):
       gp.c('set yrange[-1:1.5]')
    elif (Properties == "wpsfunction"):
       gp.c('set yrange[-0.5:1.5]')
    elif (Properties == "density" and data == "PCORE"):
       gp.c('set yrange[0:1]')
    elif (Properties == "density" and data == "CORE"):
       gp.c('set yrange[0:9]')
    elif (Properties == "potential" and data1 == "POTAE" and data2 == "POTPS"):
       gp.c('set yrange[-50: 100]')
    elif (Properties == "potential" and data1 == "POTAEC" and data2 == "POTPSC"):
       gp.c('set yrange[-500: 50]')
    else:
#       gp.c('set yrange[-50:250]')
       gp.c('set yrange[0:4]')
    gp.c(plot_command1)
    gp.c(plot_command2)
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
    Data_compare_Fig("waefunction", "WAE", "WAE_data")
    Data_compare_Fig("wpsfunction", "WPS", "WPS_data")
    Data_compare_Fig("density", "RHOV", "RHOAE_data")
    Data_compare_Fig("density", "RHOV", "RHOAE_data")
    Data_compare_Fig("density", "CORE", "CORE_data")
    Data_compare_Fig("density", "PCORE", "PCORE_data")
    Data_compare_Fig("potential", "POTAE", "POTAE_data")
    Data_compare_Fig("potential", "POTPS", "POTPS_data")
    Data_compare_Fig("potential", "POTPSC", "POTPSC_data")

#    Data_compare_Fig2("potential",  "POTAE", "POTPSC", "VASP")

if __name__ == '__main__':
    main()
    # print(__name__)

