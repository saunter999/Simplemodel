set terminal png
set output "FS.png"
set pointsize 0.01
set xtics 0.1
set ytics 0.1
set grid
plot "FS.out 0" using ($1/pi):($2/pi) ,"FS.out 1" using ($1/pi):($2/pi),"FS.out 2" using ($1/pi):($2/pi),\
"FS.out 3" using ($1/pi):($2/pi),"FS.out 4" using ($1/pi):($2/pi),"FS.out 5" using ($1/pi):($2/pi) 
