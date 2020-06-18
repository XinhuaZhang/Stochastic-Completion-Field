#!/bin/bash
GPUID="[0,1,2,3]"
Size=255
SizeRecon=128
Delta=0.25
DeltaRecon=1
Threshold=0.000001
Orientations=361
Scales=1
ThetaSigma=0.0025
ScaleSigma=0.0
MaxScale=0
Tau=10
NumR2Freq=256
PeriodR2=256
PhiFreqs=30
RhoFreqs=30
ThetaFreqs=30
ScaleFreqs=0
InitDist="[(0,-19.5,45,1),(0,20.5,135,1)]"
InitScale=1
STDR2=48
STD=15
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${NumR2Freq}_${Delta}_${Threshold}_${STDR2}.dat
HistPath=/holding_disk/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=1
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
S=1
WriteFlag=True
Iterations=10
Circle="(Circle{circleNum=8,circleRadius=30})"
Shape2D="Points(0,0)2${Circle}"
Threads=28

time stack test :IllusoryContourShapePinwheelBasis-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${WriteFlag} ${Iterations} ${Shape2D} ${Threads} +RTS -N${Threads} -RTS"

