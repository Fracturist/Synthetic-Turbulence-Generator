############################### Set Terminal ##################################

set term pngcairo enh size 1280,960 lw 3 font "Times New Roman,27"


############################# Define Variables ################################

nlw = 1.5  # lindwidth
nps = 2.5  # pointsize

purple = 1
green = 2
lblue = 3
orange = 4
yello = 5
dblue = 6
red = 7
black = 8

square = 5
circle = 7
utri = 9
dtri = 11
diamond = 13
pentagon = 15


############################## Default Reset ##################################

reset                                  # reset all

set autoscale                          # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically


############################### Plot Section ##################################

prefix = "tmp"
folder = "postProcessing/energySpectrum/0/"
file = "energySpectrum.dat"

set output prefix.".png"
set xlabel "{/symbol-Italic k} [1/m]"
set ylabel "{/Times-Italic E}({/symbol-Italic k}) [m^3/s^2]"
set title "&{aaa}"
set log xy
set xrange [10:2500]
set yrange [1e-7:1e-3]
set xtics 10
set mxtics 9
set format y "%.0e"
set key \
    bottom left at 18,0.8e-6 Left reverse spacing 1.5 samplen 4  \
    font 'Times New Roman,24'  box

plot folder.file u 1:2 w l  lc purple  lt 9  lw nlw  title "Generated"


###############################################################################
