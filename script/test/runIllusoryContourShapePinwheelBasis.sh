#!/bin/bash
GPUID="[2,3,4,5,6,7]"
Size=1024  #2047
SizeRecon=256
Delta=0.25 #0.03125
DeltaRecon=1
Threshold=0.000001
Orientations=361
Scales=1
ThetaSigma=0.0075
ScaleSigma=0.0
MaxScale=0
Tau=50
NumR2Freq=256
PeriodR2=512
PhiFreqs=30
RhoFreqs=30
ThetaFreqs=30
ScaleFreqs=0
InitDist="[(0,-19.5,45,1),(0,20.5,135,1)]"
InitScale=1 #0.149
STDR2=0.015625 #0.0078125 #0.015625 #0.015625
STDA=0
STDR=0
S=1
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}.dat
HistPath=/holding_disk/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=4
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
WriteFlag=True
Iterations=10
Circle="(Circle{circleNum=8,circleRadius=60})"
Shape2D="Points(0,0)2${Circle}"
Radius=0.75 #0.015625
Threads=28

time stack test :IllusoryContourShapePinwheelBasis-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STDA} ${STDR} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${WriteFlag} ${Iterations} ${Shape2D} ${Radius} ${Threads} +RTS -N${Threads} -RTS"

