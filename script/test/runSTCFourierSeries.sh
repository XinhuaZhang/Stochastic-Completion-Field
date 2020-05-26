#!/bin/bash
GPUID="[0,1,2,3,4,5,6,7]"
Size=64
Delta=1 #0.03125 #0.08494630713934516
Orientations=144
Scales=1
ThetaSigma=0.005
ScaleSigma=0.0
MaxScale=32
DeltaLog=0.1
Tao=25
Trails=48000
MaxTrails=4800
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=0
InitDist="[(0,0,0,1),(0,0,330,1)]"
InitScale=1
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
STD=10
STDG=1
Location="(25,4)"
Threads=28
time stack test :STCFourierSeries-test --test-arguments "${GPUID} ${Size} ${Delta} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${DeltaLog} ${Tao}  ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitDist} ${InitScale} ${HistPath} ${STD} ${STDG} ${Location} ${Threads} +RTS -N${Threads} -RTS"
