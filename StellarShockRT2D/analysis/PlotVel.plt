
#####################
# initialize
######################

set style line 1 lt 1 lw 6 lc rgb "#ff2800" # universal design red 
set style line 2 lt 1 lw 6 lc rgb "#0041ff" # universal design blue
set style line 3 lt 1 lw 6 lc rgb "#35a16B" # universal design green
set style line 4 lt 1 lw 6 lc rgb "#faf500" # universal design yellow
set style line 5 lt 1 lw 6 lc rgb "#66ccff" # universal design sky-blue,azure
set style line 6 lt 1 lw 6 lc rgb "#ff99a0" # universal design pink
set style line 7 lt 1 lw 6 lc rgb "#ff9900" # universal design orange
set style line 8 lt 1 lw 6 lc rgb "#9a0079" # universal design purple
set style line 9 lt 1 lw 6 lc rgb "#663300" # universal design brown


# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" size 640,480 
# crop 

if (exist("ifnum")==0 ) ifnum=100

ifnames = sprintf("output/velpro%05d.dat",ifnum)

command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)
flag=0
flag=system(command)
print  ifnames." found"

command = sprintf("awk 'NR==1 {print $3}' %s",ifnames)
time   = system(command)
time = time + 0.0
timeunit=" s"
timetxt = sprintf("%g",time).timeunit
print "time= ".timetxt

command = sprintf("awk 'NR==2 {print $3}' %s",ifnames)
binwidth   = system(command)
print "bin= ".binwidth
binwidth = real(binwidth) 

set xlabel "Velocity [km/s]" offset 0,0
set xrange [0:*]

ofname = sprintf("figures/veldis%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "{/Symbol D}M [M_s]" offset 0,0

bin(x) = binwidth * floor(x/binwidth)
set format y "10^{%L}"
set log y
set yrange [1e-6:*]

set style fill transparent solid 0.4 border
set boxwidth binwidth

plot ifnames u (bin($1)):2 smooth frequency with boxes lc rgb "gray" title "total"  \
,    ifnames u (bin($1)):4 smooth frequency with boxes lc rgb "blue" title "He" \
,    ifnames u (bin($1)):5 smooth frequency with boxes lc rgb "green" title "CO" \
,    ifnames u (bin($1)):6 smooth frequency with boxes lc rgb "red"  title "Fe"

print ofname." is written"


reset
set term pop
