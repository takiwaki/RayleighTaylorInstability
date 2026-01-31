##########################################
# Makefile for Blastwave Analysis
##########################################

##########################################
# Programs
##########################################

ana=Analysis.x

######################	
# complile options
######################
fc=ifort
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

# 1D radial profile
ONEPRO    = $(patsubst ${dird}/unf%.dat,${dira}/onepro%.dat,$(BIN))
DENONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/denone%.png,$(ONEPRO))
VELONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/velone%.png,$(ONEPRO))
PREONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/preone%.png,$(ONEPRO))
XCMONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/xcmone%.png,$(ONEPRO))


TWOPRO    = $(patsubst ${dird}/unf%.dat,${dira}/twopro%.dat,$(BIN))
DENTWOPNG = $(patsubst ${dira}/twopro%.dat,${dirf}/dentwo%.png,$(TWOPRO))
#VELTWOPNG = $(patsubst ${dira}/rtp%.dat,${dirf}/vel%.png,$(TWOPRO))
#PRETWOPNG = $(patsubst ${dira}/rtp%.dat,${dirf}/pre%.png,$(TWOPRO))
XCMTWOPNG = $(patsubst ${dira}/twopro%.dat,${dirf}/xcmtwo%.png,$(TWOPRO))

.PRECIOUS: ${ONEPRO}

##########################################
# Movies
##########################################

denonemovie=${dirm}/anidenone.mp4
velonemovie=${dirm}/anivelone.mp4
preonemovie=${dirm}/anipreone.mp4
xcmonemovie=${dirm}/anixcmone.mp4

dentwomovie=${dirm}/anidentwo.mp4
xcmtwomovie=${dirm}/anixcmtwo.mp4

movies = ${denonemovie} ${velonemovie} ${preonemovie} ${xcmonemovie} ${dentwomovie} ${xcmtwomovie}
##########################################
# Time evolution
##########################################

timefile=t-E.png

##########################################
# Proceadures
##########################################

all: build movie timeevolution

.PHONY: all

##########################################
# Time evolution
##########################################
.PHONY: timeevolution
timeevolution: ${timefile} ## time evolution 

${timefile}: PlotTime.plt t-prof.dat
	gnuplot PlotTime.plt

#################
# ffmpeg
#################
.PHONY: movie
movie: ${movies} ## Make movie

${denonemovie}: MakeMovie.sh ${DENONEPNG}
	./MakeMovie.sh denone

${velonemovie}: MakeMovie.sh ${VELONEPNG}
	./MakeMovie.sh velone

${preonemovie}: MakeMovie.sh ${PREONEPNG}
	./MakeMovie.sh preone

${xcmonemovie}: MakeMovie.sh ${XCMONEPNG}
	./MakeMovie.sh xcmone

${dentwomovie}: MakeMovie.sh ${DENTWOPNG}
	./MakeMovie.sh dentwo

${xcmtwomovie}: MakeMovie.sh ${XCMTWOPNG}
	./MakeMovie.sh xcmtwo

#################
# gnuplot
#################
.PHONY: 1Dsnaps
1Dsnaps: ${DENONEPNG} ${VELONEPNG} ${PREONEPNG} ${XCMONEPNG}  ## Make one-dimentional plot

${DENONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

${VELONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

${PREONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

.PHONY: 2Dsnaps
2Dsnaps: ${DENTWOPNG} ${XCMTWOPNG}  ## Make two-dimentional plot

${DENTWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh twopro Plot2D.plt

${XCMTWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh twopro Plot2D.plt


#################
# analysis
#################
run: ${ONEPRO} ${TWOPRO} ## run the analysis

count-number: ${countfile} ## count bin files to prepare the analysis

${countfile}: CountBindata.sh
	./CountBindata.sh

${ONEPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

${TWOPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

t-prof.dat: ${ana} ${BIN} ${countfile}
	./${ana}

#################
# simulation
#################

build: ${ana} ## Compile the analyis code

${ana}: Analysis.f90
	${fc} ${fopt} $< -o ${ana}

#################
# clean up
#################
.PHONY: clean allclean
clean: ## Delete intermediate files to compile
	rm -f  ${ana} *.o *.mod *~

allclean: ## Delete data
	rm -fr ${dirm} ${dirf}  ${dira} ${countfile} t-prof.dat

.PHONY: help
help: ## Display this help screen
	@grep -E '^[0-9a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'
