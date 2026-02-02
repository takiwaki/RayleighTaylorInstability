
reset

pngflag=1

if(pngflag==1) outflag=1
if(outflag==1) set terminal push

if(pngflag ==1) set terminal pngcairo enhanced font "Helvetica, 18" crop

input="t-r-pro.dat"

set pm3d map

##########################################
# Space Time Diagram
##########################################

set xlabel "Time [s]" offset 0,1
set xtic offset 0,0.5
set xtic 50
set mxtic 5
set xrange [*:*]

set ylabel "Radius [cm]" offset 0.5, 0.0
set ytic offset 0.5,0.0
set yrange [*:*]
set log y
set format y "10^{%L}"

# density
if(pngflag ==1) set output "t-r-rho.png"
if(pngflag ==1) print "output t-r-rho.png"

set palette defined (0.00 "#440154", 0.25 "#3b528b", 0.50 "#21908d", 0.75 "#5dc863", 1.00 "#fde725" ) # viridis
set title "log({/Symbol r} cm^3)"

set cbrange [*:*]

splot input  u ($1):($2):(log10($4)) notitle \

# velocity
if(pngflag ==1)set output "t-r-vr.png"
if(pngflag ==1) print "output t-r-vr.png

set palette defined (0.00 "#0d0887", 0.25 "#7e03a8", 0.50 "#cc4778", 0.75 "#f89540", 1.00 "#f0f921" ) #plasma
set title "log(v_r cm/s)"
set cbrange [*:*]
set xrange [*:*]
set yrange [*:*]

splot input  u ($1):($2):(log10($6)) notitle \


# elements
if(pngflag ==1)set output "t-r-xcm.png"
if(pngflag ==1) print "output t-r-xcm.png"

set title "X, Ni=4, CO=3, He=2, H=1"

set palette define (1.0 "black",2.0 "blue", 3.0 "green", 4.0 "red")
set cbrange [1:4]

splot input  u ($1):($2):(4*$7+3*$8+2*$9+$10) notitle \

# 1:r[cm]       2:M[Ms]       3:den[g/cm^3] 4:p[erg/cm3]  5:vel[cm/s]   6:X_Ni        7:X_CO        8:X_He        9:X_H        

##########################################
# Finalize
##########################################

reset
if(outflag==1) set terminal pop
