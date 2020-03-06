#!/bin/bash
GPUID="[7]"
Size=128
Delta=1.5 #0.03125 #0.08494630713934516
Orientations=36
Scales=1
ThetaSigma=0.04
ScaleSigma=0.0
MaxScale=32
DeltaLog=1
Tao=25
Trails=96000
MaxTrails=48000
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=0
InitDist="[(-31,0,0,1),(30,0,330,1)]"
InitScale=1
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
STD=10
Threads=24
time stack test :STCFourierSeries-test --test-arguments "${GPUID} ${Size} ${Delta} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${DeltaLog} ${Tao}  ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STD} ${Threads} +RTS -N${Threads} -RTS"
