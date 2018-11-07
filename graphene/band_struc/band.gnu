set terminal png
set output "band_kpath.png"
unset xtics
set ylabel "Energy"
set yrange [-5:5]
set label "gamma" at 0,-5.3
set label "K'" at 100,-5.3
set label "M" at 150,-5.3
set label "K" at 200,-5.3
set label "gamma" at 285,-5.3
set xzeroaxis
plot "downband.txt" u 1 w l title "downband","upperband.txt" u 1 w l title "upband"

