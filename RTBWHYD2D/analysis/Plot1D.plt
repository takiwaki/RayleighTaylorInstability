
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

unset key

# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" size 640,480 
# crop 

if (exist("ifnum")==0 ) ifnum=100

ifnames = sprintf("output/rpr%05d.dat",ifnum)

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


xnorm=1.0e10
set xlabel "Radius [10^{10} cm]" offset 0,0


##########################
# Pressure
##########################

ofname = sprintf("figures/preone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Pressure [erg/cm^3]" offset 0,0


set log 

plot  \
  ifnames u ($1/xnorm):4 w l ls 1 \

print ofname." is written"

##########################
# density
##########################_

ofname = sprintf("figures/denone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Density [g/cm^3]" offset 0,0

set log 

plot  \
  ifnames u ($1/xnorm):3 w l ls 1 \

print ofname." is written"

##########################
# velocity
##########################_

ofname = sprintf("figures/velone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "velocity [cm/s]" offset 0,0


unset log y

plot  \
  ifnames u ($1/xnorm):5 w l ls 1 \

print ofname." is written"

##########################
# X
##########################_

ofname = sprintf("figures/xcmone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Composition" offset 0,0

set key top right

unset log y

set yrange [1.0e-2:1.0]
plot  \
  ifnames u ($1/xnorm):6 w l ls 1 title "X_{Ni}"\
, ifnames u ($1/xnorm):7 w l ls 2 title "X_{CO}"\
, ifnames u ($1/xnorm):8 w l ls 3 title "X_{He}" \
, ifnames u ($1/xnorm):9 w l ls 5 title "X_{H}" \


print ofname." is written"


reset
set term pop
