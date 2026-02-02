##########################################
# Makefile for Rayleigh Taylor Instability
##########################################

##########################################
# Programs
##########################################

ana=Analysis.x

######################	
# complile options
######################
fc=ifort -diag-disable=10448
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

##########################################
# Timeseqencial data and figures
##########################################

######################
# directory
######################

dird := ../bindata
dira := output
dirf := figures
dirm := movies

countfile=control.dat

######################
# files
######################

filename = ${dird}/unf00001.dat
fileexists = $(shell ls | grep ${filename})
ifeq (${fileexists}, ${filename})
    BIN= $(shell ls ${dird}/unf*.dat)
else
    BIN= ${filename}
endif

# 2D cross section
DEN   = $(patsubst ${dird}/unf%.dat,${dira}/den%.dat,$(BIN))
DENPNG= $(patsubst ${dira}/den%.dat,${dirf}/den%.png,$(DEN))

# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

denmovie=${dirm}/aniden.mp4

2Dmovies=${denmovie}

moviefile=${2Dmovies}

##########################################
# Time evolution
##########################################

timefile=t-vsq.png t-W.png

##########################################
# Proceadures
##########################################

all: ${exe} ${moviefile} ${timefile}

.PHONY: all clean allclean

##########################################
# Time evolution
##########################################
${timefile}: TimeEvolution.plt t-prof.dat
	gnuplot TimeEvolution.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh

#################
# ffmpeg
#################

movie: ${moviefile}

${denmovie}: MakeMovie.sh ${DENPNG}
	./MakeMovie.sh den

#################
# gnuplot
#################


########################
# 2D cross section
########################

${DENPNG}: den.plt MakeSnap.sh ${DEN}
	./MakeSnap.sh den den.plt

2Dsnaps: ${DENPNG} 

#################
# analysis
#################

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

${DEN}: ${ana} ${BIN} ${countfile}
	./${ana}

#################
# simulation
#################

${ana}: Analysis.f90
	${fc} ${fopt} $< -o ${ana}

#################
# clean up
#################
clean:
	rm -f ${ana} *.o *.mod *~

allclean:
	rm -fr ${dirm} ${dirf} ${dira} ${countfile}
