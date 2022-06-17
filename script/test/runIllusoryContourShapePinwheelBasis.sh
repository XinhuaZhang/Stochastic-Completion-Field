#!/bin/bash
GPUID="[0]"
Size=512  #2047
Delta=0.25 #0.03125
SizeRecon=256
DeltaRecon=1 #0.6 # 0.16666666666666667
Threshold=0.000001
Orientations=90
Scales=1
ThetaSigma=0.018 #0.018 #0.012
ScaleSigma=0.0
Tau=9
NumTrails=1000000
BatchSize=1000000
DeltaT=0.125
PoissonWeight=0.000015
NumR2Freq=256
PeriodR2=256
DeltaFreq=0.5 #1 #0.5
PhiFreqs=32
RhoFreqs=32
ThetaFreqs=32
ScaleFreqs=32
InitDist="[(0,-19.5,45,0.5),(0,20.5,135,0.5)]"
InitScale=1 #0.149
STDR2=0.5 #0.5 #7.295125 #0.0078125 #0.015625 #0.015625
STDTheta=0.3 #7.295125
STDR=0.5 #1.5
S=1
HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${STDR2}.dat
HistFileNameCorner=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${STDR2}_${PoissonWeight}.dat
# HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${DeltaT}.dat
HistPath=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileName}
HistPathCorner=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileNameCorner}
NumBatchR2=16
NumBatchR2Freqs=1
NumBatchOri=1
WriteFlag=True
Iterations=20
Circle="(Circle{circleNum=8,circleRadius=$1})"
KoffkaCross1="(KoffkaCross{koffkaCrossW=24,koffkaCrossD=48})"
KoffkaCross2="(KoffkaCross{koffkaCrossW=12,koffkaCrossD=48})"
Shape2D="[Points(0,0)2${KoffkaCross1},Points(0,0)2${KoffkaCross2}]"
Shape2D="[Points(0,0)2${Circle}]"
Radius=1.5 #.015625
Threads=20

time stack test :IllusoryContourShapePinwheelBasis-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${Tau} ${NumTrails} ${DeltaT} ${PoissonWeight} ${NumR2Freq} ${PeriodR2} ${DeltaFreq} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${HistPathCorner} ${STDR2} ${STDTheta} ${STDR} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${WriteFlag} ${Iterations} ${Shape2D} ${Radius} ${Threads} +RTS -N${Threads} -RTS"

