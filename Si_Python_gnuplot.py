#!/home/jun_jiang/.local/bin/python
#!coding = utf-8
import os
import time

# g = PyGnuplot.Gnuplot()   

#import PyGnuplot as gp 
from PyGnuplot import gp 
import numpy as np 

# X=np.crange(10)
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
LDA-PW scalarrelativistic loggridv4 5001 9.2 5.20 loggderivrange -25. 25. 2500
3 3 0 0 0 0
3 1 2
0 0 0
c
c
v
c
v
1
1.90  1.41   1.44    1.505 # 1.85 1.51 1.54 1.51
y
-0.35   # -0.70 能量往大调压制波函数
n
y
0.50    # 0.05 能量往大调压制波函数
n
vasprrkj  vanderbiltortho  Besselshape  #Vanderbilt vasprrkj  vanderbiltortho  Besselshape
1 3   bessel  #1 3   MTROULLIER
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

    plot_command1 = '"'+type+'_'+data +'" '+'u 1:2 w p pt 3 ps 1.2 title "'+TYPE +'", '
    plot_command2 = '"'+type+'_'+data +'" '+'u 1:3 w l lt 5 lw 1.2 title "'+TYPE +'\n"'
#    plot_command_ATOM = "replot "+'"ATOM_'+data +'" '+'u 1:2 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + data +'"'
#    Fig_name = type+"_"+name+".eps"
    Fig_name = type+"_"+name+".ps"

    Fig = gp()
    Fig.c(Fig_title)
    if (Properties == "wavefunction"):
       Fig.c('set xrange[0:1.5]')
       Fig.c('set yrange[-1.5:1.5]')

    Fig.c("plot "+ plot_command1 + plot_command2)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
#    Fig.ps(Fig_name)
    Fig.ps(Fig_name)

def Data_compare_Fig(Properties, data, name):
    TYPE1 = "VASP-"+data
    TYPE2 = "ATOM-"+data

    plot_command_VASP = '"VASP_'+data +'" '+'u 1:2 w p pt 1 lw 1.0 title "'+TYPE1+'", '
    plot_command_ATOM = '"ATOM_'+data +'" '+'u 1:2 w l lt 3 lw 1.5 title "'+TYPE2+'"'
#    plot_command_ATOM = "replot "+'"VASP_'+data +'" '+'u 1:3 w l lt 5 lw 1.8 title "'+TYPE2+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + name +'"'
#    Fig_name = name+".eps"
    Fig_name = name+".ps"

    Fig = gp()
    Fig.c(Fig_title)
    Fig.c('set xrange[0:3.0]')
    if (Properties == "waefunction"):
       Fig.c('set yrange[-1:1.5]')
    elif (Properties == "wpsfunction"):
       Fig.c('set yrange[-0.5:1.5]')
    elif (Properties == "density" and data == "PCORE"):
       Fig.c('set yrange[0:1]')
    elif (Properties == "density" and data == "CORE"):
       Fig.c('set yrange[0:9]')
    elif (Properties == "potential" and data == "POTAE"):
       Fig.c('set yrange[-50: 5000]')
    elif (Properties == "potential" and data == "POTPS"):
       Fig.c('set yrange[-50:20]')
    elif (Properties == "potential" and data == "POTPSC"):
       Fig.c('set yrange[-500:100]')
    else:
#       gp.c('set yrange[-50:250]')
       Fig.c('set yrange[0:4]')
    Fig.c("plot" + plot_command_VASP + plot_command_ATOM)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
    Fig.ps(Fig_name)
#    Fig.pdf(Fig_name)

def Data_compare_Fig2(Properties, data, name):
    TYPE1 = "VASP-"+data
    TYPE2 = "ATOM-"+data
    TYPE3 = "VASP-"+"WAV-CAL"

    plot_command1 = '"VASP_'+data +'" '+'u 1:2 w p pt 3 lw 1.2 title "'+TYPE1+'", '
    plot_command2 = '"ATOM_'+data +'" '+'u 1:2 w l lt 5 lw 1.8 title "'+TYPE2+'", '
    plot_command3 = '"VASP_'+data +'" '+'u 1:3 w l lt 7 lw 1.8 title "'+TYPE3+'"'
#    print(plot_command_VASP)
#    exit ()
    Fig_title = "set title "+ '"' + name +'"'
#    Fig_name = name+".eps"
    Fig_name = name+".ps"
    Fig = gp()

    Fig.c(Fig_title)
    Fig.c('set xrange[0:3.0]')
    if (Properties == "waefunction"):
       Fig.c('set yrange[-1:1.5]')
    elif (Properties == "wpsfunction"):
       Fig.c('set yrange[-0.5:1.5]')
    elif (Properties == "density" and data == "PCORE"):
       Fig.c('set yrange[0:1]')
    elif (Properties == "density" and data == "CORE"):
       Fig.c('set yrange[0:9]')
    elif (Properties == "potential" and data == "POTAE"):
       Fig.c('set yrange[-50: 100]')
    elif (Properties == "potential" and data == "POTPS"):
       Fig.c('set yrange[-50: 20]')
#    elif (Properties == "potential" and data1 == "POTAE" and data2 == "POTPS"):
#       gp.c('set yrange[-50: 100]')
#    elif (Properties == "potential" and data1 == "POTAEC" and data2 == "POTPSC"):
#       gp.c('set yrange[-500: 50]')
    else:
#       gp.c('set yrange[-50:250]')
       Fig.c('set yrange[0:4]')
    Fig.c("plot" + plot_command1+plot_command2+plot_command3)
    # gp.c('plot "VASP_WAE" u 1:2 w p pt 3 ps 2 title "VASP-AE"')
    # gp.c(plot_command_ATOM)
    # gp.p('myfigure.ps')
    Fig.ps(Fig_name)
#    Fig.pdf(Fig_name)

def main():
    Write_Atom_in('Si', 14)

    os.system("~/Softs/atompaw_test/src/atompaw < in > out")
    #exit ()
#    time.sleep(5)


#    Data_Fig("potential", "VASP", "POTAE", "POTAE-data")
    Data_compare_Fig("waefunction", "WAE", "WAE-data")
    Data_compare_Fig("wpsfunction", "WPS", "WPS-data")
#    Data_compare_Fig("density", "RHOV", "RHOAER-ata")
#    Data_compare_Fig("density", "RHOV", "RHOAER-ata")
    Data_compare_Fig("density", "CORE", "CORE-data")
    Data_compare_Fig("density", "PCORE", "PCORE-data")
    Data_compare_Fig2("potential", "POTAE", "POTAE-data")
    Data_compare_Fig2("potential", "POTPS", "POTPS-data")
    Data_compare_Fig("potential", "POTPSC", "POTPSC-data")

#    Data_compare_Fig2("potential",  "POTAE", "POTPSC", "VASP")
    os.system("evince *-data.ps")
    exit()

if __name__ == '__main__':
    main()
    # print(__name__)

