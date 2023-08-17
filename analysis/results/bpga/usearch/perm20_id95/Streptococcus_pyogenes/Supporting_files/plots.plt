#line 4303
set style data dots
set style function lines
set xzeroaxis lt 6 lc rgb "grey" lw 3
set terminal unknown
set xrange [0.95:]

set terminal postscript enhanced size 10,7
set offset graph 0.0, graph 0.1, graph 0.1, graph 0.1
set output '| ps2pdf - Core_Pan_Plot.pdf'
set title "Pan and Core Genome Plot" font "arial bold, 18" textcolor rgb"blue"
set xlabel "Number of Genomes\n\n" font "arial bold, 16" textcolor rgb"black"
set ylabel "\nNumber of Gene Families\n" font "arial bold, 16" textcolor rgb"black"
set key outside horizontal bottom title "" font "arial , 14" textcolor rgb"black"

set bars 3.0
set boxwidth  -2
set style fill empty
set xtics 1 font "arial , 14" textcolor rgb"black"
set ytics 500 font "arial , 14" textcolor rgb"black"

plot 'core_pan_box.txt' using 6:2:1:5:4 with candlesticks lc rgb"black" title 'Pan \& Core Genome' whiskerbars , \
''         using 6:3:3:3:3 with candlesticks lt 6 lc rgb"red" title 'Median Values'

unset output

exit gnuplot;

