#!/bin/bash
GPUID="[0,1,2,3,4,5,6,7]"
Size=1024
SizeRecon=256
Delta=0.25
DeltaRecon=1
Threshold=0.00001
Orientations=360
Scales=1
ThetaSigma=0.2 #0.0005
ScaleSigma=0.0
Tau=7.5 #9.5
NumR2Freq=256
PeriodR2=256
EnvelopePeriod=363
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=20
InitDist="[(0,0,0,0.5),(0,0,0,8)]"
InitScale=1 #0.149
STDR2=1
STD=15
S=1
DeltaT=0.5
HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${DeltaT}.dat
# HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${NumR2Freq}_${Delta}_${Threshold}_${STDR2}.dat
HistPath=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=32
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
Radius=8
Threads=20

time stack test :STCPinwheel-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${Radius} ${DeltaT} ${Threads} +RTS -N${Threads} -RTS"

