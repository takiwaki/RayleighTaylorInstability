
##########################################
# Makefile for Turbulence Studies
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
ONEPRO    = $(patsubst ${dird}/unf%.dat,${dira}/rpr%.dat,$(BIN))
DENONEPNG = $(patsubst ${dira}/rpr%.dat,${dirf}/den%.png,$(ONEPRO))
VELONEPNG = $(patsubst ${dira}/rpr%.dat,${dirf}/vel%.png,$(ONEPRO))
PREONEPNG = $(patsubst ${dira}/rpr%.dat,${dirf}/pre%.png,$(ONEPRO))
XCMONEPNG = $(patsubst ${dira}/rpr%.dat,${dirf}/xcm%.png,$(ONEPRO))


TWOPRO    = $(patsubst ${dird}/unf%.dat,${dira}/rtp%.dat,$(BIN))
DENTWOPNG = $(patsubst ${dira}/rtp%.dat,${dirf}/dnt%.png,$(TWOPRO))
#VELTWOPNG = $(patsubst ${dira}/rtp%.dat,${dirf}/vel%.png,$(TWOPRO))
#PRETWOPNG = $(patsubst ${dira}/rtp%.dat,${dirf}/pre%.png,$(TWOPRO))
XCMTWOPNG = $(patsubst ${dira}/rtp%.dat,${dirf}/xct%.png,$(TWOPRO))


# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

denonemovie=${dirm}/aniden_one.mp4
velonemovie=${dirm}/anivel_one.mp4
preonemovie=${dirm}/anipre_one.mp4
xcmonemovie=${dirm}/anixcm_one.mp4

dentwomovie=${dirm}/aniden_two.mp4
xcmtwomovie=${dirm}/anixcm_two.mp4

movies = ${denonemovie} ${velonemovie} ${preonemovie} ${xcmonemovie} ${dentwomovie} ${xcmtwomovie}
##########################################
# Time evolution
##########################################

timefile=t-E.png

##########################################
# Proceadures
##########################################

all: build movie timeevolution

.PHONY: all clean allclean

##########################################
# Time evolution
##########################################
.PHONY: timeevolution
timeevolution: ${timefile} ## time evolution 

${timefile}: t-E.plt t-prof.dat
	gnuplot t-E.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh

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
.PHONY: onedim
onedim: ${DENONEPNG} ${VELONEPNG} ${PREONEPNG} ${XCMONEPNG}  ## Make one-dimentional plot

${DENONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh rpr Plot1D.plt

${VELONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh rpr Plot1D.plt

${PREONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh rpr Plot1D.plt

.PHONY: twodim
twodim: ${DENTWOPNG} ${XCMTWOPNG}  ## Make two-dimentional plot

${DENTWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh rtp Plot2D.plt

${XCMTWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh rtp Plot2D.plt


#################
# analysis
#################
run: ${ONEPRO} ${TWOPRO} ## run the analysis

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

${ONEPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

${TWOPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

${TOT}: ${ana} ${BIN} ${countfile}
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
clean: ## Delete intermediate files to compile
	rm -f  ${ana} *.o *.mod *~

allclean: ## Delete data
	rm -fr ${dirm} ${dirf}  ${dira} ${countfile}

.PHONY: help
help: ## Display this help screen
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'
