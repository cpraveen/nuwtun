
set xlabel 'Number of iterations'
set ylabel 'C_l'
set y2label 'C_d'
set ytics nomirror
set y2tics
set term postscript enhanced
set out 'clcd.eps'
p 'fort.18' u 1:3 t 'C_l' w l lw 2 axis x1y1, \
  'fort.18' u 1:4 t 'C_d' w l lw 2 axis x1y2

set term x11
set out

unset y2label
unset y2tics

set logscale y
set xlabel 'Number of iterations'
set ylabel 'Residue'
set term postscript enhanced
set out 'res.eps'
p 'fort.17' u 1:3 t 'Density' w l lw 2, \
  'fort.17' u 1:4 t 'x-momen' w l lw 2, \
  'fort.17' u 1:5 t 'y-momen' w l lw 2, \
  'fort.17' u 1:6 t 'z-momen' w l lw 2, \
  'fort.17' u 1:7 t 'Energy ' w l lw 2

set term x11
set out
