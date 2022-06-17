#!/bin/bash
GPUID="[0]"
Size=512  #2047
Delta=0.25 #0.03125
SizeRecon=256
DeltaRecon=1 # 0.16666666666666667
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
DeltaFreq=0.5
PhiFreqs=32
RhoFreqs=32
ThetaFreqs=32
ScaleFreqs=32
InitScale=1 #0.149
STDR2=0.5 #0.5 #7.295125 #0.0078125 #0.015625 #0.015625
STDTheta=0.5 #0.3 #7.295125
STDR=0.5 #0.5 #1.5
S=1
HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${STDR2}.dat
HistFileNameCorner=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${STDR2}_${PoissonWeight}.dat
# HistFileName=${ThetaSigma}_${ScaleSigma}_${Tau}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}_${Size}_${Delta}_${PeriodR2}_${Threshold}_${S}_${DeltaT}.dat
HistPath=/scratch/users/xinhuazhang/Pinwheel/Z2S1R+/${HistFileName}
NumBatchR2=8
NumBatchR2Freqs=1
NumBatchOri=1
WriteFlag=True
Iterations=200
Radius=1.25 #.015625
EdgeFilePath=/home/xinhuazhang/Downloads/avocados1.edge
EdgeFilePath1=/home/xinhuazhang/Downloads/Contours/N170322-001-999_9.txt
EdgeFilePath2=/home/xinhuazhang/Downloads/Contours/N180522-002-999_4.txt
EdgeFilePath3=/home/xinhuazhang/Downloads/Contours/N170322-001-999_8.txt
# EdgeFilePath1=/home/xinhuazhang/Downloads/x_vals.txt

# ScaleFactor=0.6 #0.75
ScaleFactor=0.75
NumNoisePoint=20
NoisePointRange="(-100,100)"
Threads=16
export ACCELERATE_LLVM_NATIVE_THREADS=${Threads}

time stack test :IllusoryContourEdgeCorner-test --test-arguments "${GPUID} ${Size} ${Delta} ${Threshold} ${SizeRecon} ${DeltaRecon} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${Tau} ${NumTrails} ${DeltaT} ${PoissonWeight} ${NumR2Freq} ${PeriodR2} ${DeltaFreq} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitScale} ${HistPath} ${STDR2} ${STDTheta} ${STDR} ${NumBatchR2} ${NumBatchR2Freqs} ${NumBatchOri} ${BatchSize} ${S} ${WriteFlag} ${Iterations} ${Radius} ${EdgeFilePath1} ${EdgeFilePath2} ${EdgeFilePath3} ${ScaleFactor} ${NumNoisePoint} ${NoisePointRange} ${Threads} +RTS -N${Threads} -RTS"

