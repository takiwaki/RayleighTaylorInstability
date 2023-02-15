#####################
# initialize
######################

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set view map
unset key
set size ratio -1

##########################################
# parameters
##########################################

srange=0.5

##########################################
# mainloop
##########################################

# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" size 300,500
# crop 

if (exist("ifnum")==0 ) ifnum=50

ifnames = sprintf("output/den%05d.dat",ifnum)

command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)
flag=0
flag=system(command)
print  ifnames." found"

command = sprintf(" head -n 1 %s | sed 's/#//' ",ifnames)
time   = system(command)
print "time= ".time

# Position of color bar
set colorbox horizontal user origin 0.3, 0.87 size 0.4, 0.04
#set cbtics 1.0
set cbtics offset 0,3.2

set size 1.0

set origin 0.05,0.0
set xlabel "X" offset 0,0
set xtics 0.25
set ylabel "Y" offset 0,0
set ytics 0.25

set xrange [-0.25:0.25]
set yrange [-0.75:0.75]

##########################
# Density
##########################

#set title "Density"

#maxcb=20.0
#set cbrange [-maxcb:maxcb]
set cbtics 1.0
set cbrange [*:*]
set cbrange [0.9:2.1]

ofname = sprintf("figures/den%05d.png",ifnum)
set output ofname


set label 1 time at screen 0.45, screen 0.845

set palette define (-1.0 "blue", 0.0 "white", 1.0 "red")

splot  \
  ifnames u ($2):($1):3 w pm3d  \

unset label 1

reset
set term pop
