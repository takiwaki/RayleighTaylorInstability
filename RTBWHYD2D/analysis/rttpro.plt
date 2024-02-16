



####################
# Output control
####################

# Format of the output
pngflag=1

# OUTPUT PNG
if (pngflag==1) set term push
if (pngflag==1) set term pngcairo enhanced font "Helvetica, 12" size 550,500


##########################################
# initialize
##########################################

# color palette
#set palette rgbformulae 21,22,23 # black red yellow white
#set palette functions sqrt(gray), gray**3, sin(gray*2*pi) # black blue red yellow
#set palette rgb 33,13,10 # blue-green-yellow-red
#set palette define (0.0 "black",0.05 "blue", 0.5 "red", 0.75 "yellow", 1.0 "yellow")
#set palette define (-1.0 "blue", 0.0 "white", 1.0 "red")

# Lines
set style line 10 lt  1 lw 4 lc rgb "white" 

####################
# Input control
####################

# ifnum : Input File NUMber
if (exist("ifnum")==0 ) ifnum=100

print ifnum
# Scalar
ifnames = sprintf("output/rtp%05d.dat",ifnum)

# Stop Ploting
command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)

flag=0
flag=system(command)
if(flag ne "0") print ifnames." not found"; quit

# Extract Time
print  ifnames." found"

command = sprintf("awk 'NR==1 {print $3}' %s",ifnames)
time   = system(command)
print time
time = time + 0.0
timeunit=" s"
timetxt = sprintf("%g",time).timeunit
print "time=".timetxt

# Showing Time
set label timetxt at screen 0.65, screen 0.85


##########################################
# parameters
##########################################

# Range of the plot [Rs]
#stats ifnames u 1:2
#srange=STATS_max_x

command =  sprintf("awk 'NR==2 {print $3}' %s",ifnames)
rmax = system(command)
rmax = rmax +0.0
print sprintf("%g",rmax)." [cm]"
xnorm=1.0e10
set xlabel "X [10^{10} cm]" offset 0,1.0
set ylabel "Z [10^{10} cm]" offset 1.0,0
srange = rmax/xnorm


####################
# Annotation
####################

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set size ratio -1
set view map
unset key

set size 1.0, 1.0
unset origin

# vertical and horizontal axis
set origin 0.0,0.0
set xtics offset 0,0.7
#set xtics 50

set ytics offset 0,0.7
#set ytics 50

####################
# Plot
####################
ofname = sprintf("figures/dnt%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname

# Position of color bar
set colorbox horizontal user origin 0.235, 0.87 size 0.5, 0.04
set cbtics offset 0,3.2

# Range of color bar
# range of the variable
#cmin=0
#cmax=5
#set cbrange [cmin:cmax]

# Main plot
splot [-srange:srange][-srange:srange] \
  ifnames u ( $1/xnorm*sin($2)):($1/xnorm*cos($2)):($1/xnorm<srange?($3):NaN) w pm3d \
, ifnames u (-$1/xnorm*sin($2)):($1/xnorm*cos($2)):($1/xnorm<srange?($3):NaN) w pm3d \

unset label


####################
# Plot
####################
ofname = sprintf("figures/xct%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname

set palette define (1.0 "black",2.0 "blue", 3.0 "green", 4.0 "red")
set cbrange [1:4]

splot [-srange:srange][-srange:srange] \
  ifnames u ( $1/xnorm*sin($2)):($1/xnorm*cos($2)):($1/xnorm<srange?(4*$6+3*$7+2*$8+$9):NaN) w pm3d \
, ifnames u (-$1/xnorm*sin($2)):($1/xnorm*cos($2)):($1/xnorm<srange?(4*$6+3*$7+2*$8+$9):NaN) w pm3d \

unset label


##########################################
# finalize
##########################################

unset output
if(pngflag==1)set term pop

undefine ifnum