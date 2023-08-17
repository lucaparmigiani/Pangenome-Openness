set fit quiet
set fit logfile 'plots_default.log'
set boxwidth  0.3
set key outside horizontal bottom right title "" font "arial bold, 16" textcolor rgb"black"
set xtics 5 font "arial bold, 16" textcolor rgb"black"
set ytics 5000 font "arial bold, 16" textcolor rgb"black"
set offset graph 0.0, graph 0.0, graph 0.1, graph 0.0
set xrange [0:]
set style fill solid border -1

set terminal postscript enhanced size 10,7
set output "| ps2pdf - Histogram.pdf"
set title "\nDistribution of Gene Families\n" font "arial bold, 18" textcolor rgb"blue"
set xlabel "\nNumber of Genomes\n" font "arial bold, 16" textcolor rgb"black"
set ylabel "\nNumber of Gene Families\n" font "arial bold, 16" textcolor rgb"black"
plot "histogram.txt" u 1:2 w boxes lc rgb"black" notitle
unset output

set terminal postscript enhanced size 10,7
set ytics 500 font "arial bold, 16" textcolor rgb"black"
set offset graph 0.0, graph 0.0, graph 0.1, graph 0.0
set xrange [1:]

set output "| ps2pdf - New_Genes_Plot.pdf" 
set title "\nNumber of New Genes\n" font "arial bold, 18" textcolor rgb "blue"
set xlabel "\nGenome\n" font "arial bold, 16" textcolor rgb "black"
set ylabel "\nNumber of New Genes\n" font "arial bold, 16" textcolor rgb "black"
plot "new_genes_count.txt" u 1:2 w boxes lc rgb"black" notitle
unset output
set offset graph 0.0, graph 0.1, graph 0.1, graph 0.1
set style data dots

set xrange [0.95:]
set style function lines 
set pointsize 2
set style fill empty border -1
set xzeroaxis lt 6 lc rgb "grey" lw 3
set terminal unknown
f(x) = a*x**b
b = 1
a = 10000
fit f(x) 'pan_default.txt' using 1:2 via a, b
plot 'pan_default.txt' using 1:2 with dots notitle ,\
f(x) lc rgb "orange" lw 5 title '                                               Pan genome'  

f1(x) = c*exp(d*x)
d = -1
c = 10000
fit f1(x) 'core_default.txt' using 1:2 via c, d
replot 'core_default.txt' using 1:2 with dots notitle ,\
f1(x) lc rgb "purple" lw 5 title '                                               Core genome'  

replot "pan_default.txt" using 1:2 with points ps 1.5 pt 9  title '                                               Total gene families'  

unset output
set terminal postscript enhanced size 10,7
set title "\nCore-Pan Plot\n" font "arial bold, 20" textcolor rgb"blue"
set xlabel "\nNumber of Genomes\n" font "arial bold, 18" textcolor rgb"black"
set ylabel "\nNumber of Gene Families\n" font "arial bold, 18" textcolor rgb"black"
set xtics 5 font "arial bold, 16" textcolor rgb"black"
set ytics 5000 font "arial bold, 16" textcolor rgb"black"
set output '| ps2pdf - Default_Core_Pan_Plot.pdf'

replot "core_default.txt" using 1:2 with points ps 1.5 pt 9 lc rgb"pink" title '                                               Core gene families'
unset output
exit gnuplot;

