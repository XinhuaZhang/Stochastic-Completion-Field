#!/bin/bash
GPUID="[0,1,2,3,4,5,6,7]"
Size=1024  #2047
SizeRecon=128
Delta=0.25 #0.03125
DeltaRecon=1
Threshold=0.000001
Orientations=361
Scales=1
ThetaSigma=0.15
ScaleSigma=0.0
Tau=10
NumTrails=1000000
BatchSize=1000000
DeltaT=0.25
NumR2Freq=128
PeriodR2=128
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=20
InitDist="[(0,-19.5,45,0.5),(0,20.5,135,0.5)]"
InitScale=1 #0.149
STDR2=1 #7.295125 #0.0078125 #0.015625 #0.015625
STDTheta=0.5 #7.295125
STDR=0.5
S=1
HistFileName=${ThetaSigma}_${ScaleSigma}_0_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}.dat
# HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${DeltaT}.dat
HistPath=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=4
NumBatchR2Freqs=1
NumBatchOri=1
WriteFlag=True
Iterations=10
Circle="(Circle{circleNum=8,circleRadius=20})"
Shape2D="Points(0,0)2${Circle}"
Radius=1.5 #.015625
Threads=14

time stack test :IllusoryContourShapePinwheelBasis-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${Tau} ${NumTrails} ${DeltaT} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STDTheta} ${STDR} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${WriteFlag} ${Iterations} ${Shape2D} ${Radius} ${Threads} +RTS -N${Threads} -RTS"

