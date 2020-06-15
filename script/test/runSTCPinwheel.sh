#!/bin/bash
GPUID="[4,5,6,7]"
Size=257
SizeRecon=64
Delta=0.125
DeltaRecon=0.5
Threshold=0.00001
Orientations=37
Scales=1
ThetaSigma=0.0005
ScaleSigma=0.0
MaxScale=0
Tau=9.5
NumR2Freq=128
PeriodR2=64
PhiFreqs=50
RhoFreqs=50
ThetaFreqs=50
ScaleFreqs=0
InitDist="[(0,0,0,1)]"
InitScale=0.149
STDR2=48
STD=20
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${NumR2Freq}_${Delta}_${Threshold}_${STDR2}.dat
HistPath=/holding_disk/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=1
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
S=1
Threads=28

time stack test :STCPinwheel-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${Threads} +RTS -N${Threads} -RTS"

