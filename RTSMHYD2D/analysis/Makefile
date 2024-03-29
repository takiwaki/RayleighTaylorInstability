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
fc=ifort -extend-source
#fopt=-g -traceback -O2
fopt=-g -traceback -check all -fpe0

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

timefile=t-Gamma.png

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

2Dsnaps: ${VORPNG} ${KINPNG} 

${DENPNG}: den.plt MakeSnap.sh ${DEN}
	./MakeSnap.sh den den.plt

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
