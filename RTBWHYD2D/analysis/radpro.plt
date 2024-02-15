#####################
# initialize
######################

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

command = sprintf(" head -n 1 %s | sed 's/#  time_s= *//' ",ifnames)
time   = system(command)
timeunit=" s"
print "time= ".time.timeunit


set xlabel "Radius [R_s]" offset 0,0


##########################
# Pressure
##########################

ofname = sprintf("figures/pre%05d.png",ifnum)
set output ofname

timetxt = time.timeunit
set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Pressure [erg/cm^3]" offset 0,0

plot  \
  ifnames u ($1):3 w l lw 6 \

print ofname." is written"

##########################
# density
##########################_

ofname = sprintf("figures/den%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Density [g/cm^3]" offset 0,0

plot  \
  ifnames u ($1):2 w l lw 6 \

print ofname." is written"

##########################
# velocity
##########################_

ofname = sprintf("figures/vel%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "velocity [km/s]" offset 0,0

plot  \
  ifnames u ($1):4 w l lw 6 \

print ofname." is written"

##########################
# X
##########################_

ofname = sprintf("figures/xcm%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Composition" offset 0,0

set key top right

set yrange [0.0:1.0]
plot  \
  ifnames u ($1):5 w l lw 6 title "X_{Ni}"\
, ifnames u ($1):6 w l lw 6 title "X_{core}"\
, ifnames u ($1):7 w l lw 6 title "X_{H}" \


print ofname." is written"


reset
set term pop
