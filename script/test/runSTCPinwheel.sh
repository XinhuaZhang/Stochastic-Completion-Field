#!/bin/bash
GPUID="[4,5,6,7]"
Size=512
SizeRecon=128
Delta=0.25
DeltaRecon=1
Threshold=0.0000001
Orientations=361
Scales=1
ThetaSigma=0.005 #0.0005
ScaleSigma=0.0
MaxScale=0
Tau=20 #9.5
NumR2Freq=256
PeriodR2=256
EnvelopePeriod=363
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=0
InitDist="[(0,0,0,1),(0,0,90,1)]"
InitScale=1 #0.149
STDR2=1
STD=15
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${NumR2Freq}_${Delta}_${Threshold}_${STDR2}.dat
HistPath=/holding_disk/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=4
NumBatchR2Freqs=1
NumBatchOri=1
BatchSize=32
S=1
Threads=28

time stack test :STCPinwheel-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tau} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${EnvelopePeriod} ${Threads} +RTS -N${Threads} -RTS"

