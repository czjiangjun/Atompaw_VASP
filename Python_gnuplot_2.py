from subprocess import *  
gnuplot = Popen('gnuplot',stdin = PIPE, stderr=PIPE, stdout=PIPE)  
#gnuplot.stdin .write(b'set terminal jpeg \n')  
#gnuplot.stdin .write(b'set output &#39;plot.jpg&#39; \n')  #存到图片时候不会跳出gnuplot界面
gnuplot.stdin .write(b"set xlabel 'sample' \n")  
gnuplot.stdin .write(b"set ylabel 'value'\n")  
gnuplot.stdin .write(b"set title 'svm' \n")  
gnuplot.stdin .write(b"plot 'VASP_POTPS' using 1:2 with p pointsize 0.5  linetype 8,")  
gnuplot.stdin .write(b" 'VASP_POTPS' using 1:3 with p pointsize 0.5  linetype 7 \n")  
gnuplot.stdin .flush()

input('Press the Return key to quit: ') #暂停一下, 防止gnuplot界面自动关闭
gnuplot.stdin .close() #这条语句执后gnuplot界面关闭
