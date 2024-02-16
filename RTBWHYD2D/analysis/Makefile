
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
fc=ifort -extend-source
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
RPR   = $(patsubst ${dird}/unf%.dat,${dira}/rpr%.dat,$(BIN))
DENPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/den%.png,$(RPR))
VELPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/vel%.png,$(RPR))
PREPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/pre%.png,$(RPR))
XCMPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/xcm%.png,$(RPR))


RTP    = $(patsubst ${dird}/unf%.dat,${dira}/rtp%.dat,$(BIN))
DENTPNG= $(patsubst ${dira}/rtp%.dat,${dirf}/dnt%.png,$(RTP))
#VELTPNG= $(patsubst ${dira}/rtp%.dat,${dirf}/vel%.png,$(RTP))
#PRETPNG= $(patsubst ${dira}/rtp%.dat,${dirf}/pre%.png,$(RTP))
XCTTPNG= $(patsubst ${dira}/rtp%.dat,${dirf}/xct%.png,$(RTP))


# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

denmovie=${dirm}/aniden.mp4
velmovie=${dirm}/anivel.mp4
premovie=${dirm}/anipre.mp4
xcmmovie=${dirm}/anixcm.mp4


dntmovie=${dirm}/anidnt.mp4
xctmovie=${dirm}/anixct.mp4

movies = ${denmovie} ${velmovie} ${premovie} ${dntmovie}
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

${denmovie}: MakeMovie.sh ${DENPNG}
	./MakeMovie.sh den

${velmovie}: MakeMovie.sh ${VELPNG}
	./MakeMovie.sh vel

${premovie}: MakeMovie.sh ${PREPNG}
	./MakeMovie.sh pre

${xcmmovie}: MakeMovie.sh ${XCMPNG}
	./MakeMovie.sh xcm

${dntmovie}: MakeMovie.sh ${DENTPNG}
	./MakeMovie.sh dnt

${xctmovie}: MakeMovie.sh ${XCTTPNG}
	./MakeMovie.sh xct

#################
# gnuplot
#################

${DENPNG}: MakeSnap.sh radpro.plt ${RPR}
	./MakeSnap.sh rpr radpro.plt

${velPNG}: MakeSnap.sh radpro.plt ${RPR}
	./MakeSnap.sh rpr radpro.plt

${PREPNG}: MakeSnap.sh radpro.plt ${RPR}
	./MakeSnap.sh rpr radpro.plt

${DENTPNG}: MakeSnap.sh rttpro.plt ${RPR}
	./MakeSnap.sh rtp rttpro.plt

${XCMTPNG}: MakeSnap.sh rttpro.plt ${RPR}
	./MakeSnap.sh rtp rttpro.plt



#################
# analysis
#################
run: ${RPR} ## run the analysis

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

${RPR}: ${ana} ${BIN} ${countfile}
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
