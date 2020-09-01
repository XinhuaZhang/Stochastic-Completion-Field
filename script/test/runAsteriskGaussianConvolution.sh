#!/bin/bash
GPUID="[0]"
Size=256
SizeRecon=256
Delta=0.25
DeltaRecon=0.5
Threshold=0.000001
Orientations=85
Scales=1
ThetaSigma=0.02 #0.0005
ScaleSigma=0.0
MaxScale=0
Tau=20 #9.5
NumR2Freq=256
PeriodR2=64
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=20
InitDist="[(0,0.5,0,1),(0,20,150,2)]"
InitScale=1 #0.149
STDR2=1
STD=15
S=1
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}.dat
HistPath=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=4
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
Radius=8
Threads=40

time stack test :AsteriskGaussianConvolution-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${Radius} ${Threads} +RTS -N${Threads} -RTS"

