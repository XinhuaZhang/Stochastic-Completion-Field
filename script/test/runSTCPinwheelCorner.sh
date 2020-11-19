#!/bin/bash
GPUID="[0,1]"
Size=512
SizeRecon=256
Delta=0.25
DeltaRecon=1
Threshold=0.00000001
Orientations=90
Scales=1
ThetaSigma=0.015
ScaleSigma=0.0
Tau=10 #9.5
NumR2Freq=256
PeriodR2=256
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=0
InitDist="[(0,-8,0,1),(-8,0,90,1)]"
InitScale=1 #0.149
STDR2=0.5
STD=15
S=1
DeltaT=0.5
Weight=0.00001
Pos="(0,0)"
# HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${DeltaT}.dat
HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${NumR2Freq}_${Delta}_${Threshold}_${STDR2}_${Weight}.dat
HistPath=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=32
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
Radius=8
Threads=10

time stack test :STCPinwheelCorner-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${Radius} ${DeltaT} ${Weight} ${Pos} ${Threads} +RTS -N${Threads} -RTS"

