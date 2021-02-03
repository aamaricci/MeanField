reset
set terminal postscript eps enhanced color "Times-Roman" 20
set output "|ps2pdf -DEPSCrop - mf_diagram.pdf"
set yrange [0:]      #<< adjust the Yrange
set key samplen 2
set grid
unset key
unset colorbox
set xlabel '{/Symbol l}_{SO}/t'
set ylabel 'U/t'
plot 'params.run' u 2:1:8 with points palette pointsize 1 pointtype 5

