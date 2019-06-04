#!/bin/bash
Size=128
Orientations=16
Scales=8
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.000001
TaoDecay=100
Trails=1000000
MaxTrails=100000
Theta0Freqs=10
ThetaFreqs=10
Scale0Freqs=0
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${TaoDecay}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0Amodal/${HistFileName}
NumInteration=5
WriteSourceFlag=True
Sigma=0.05
InitDist="[R2S1RPPoint(15,10,330,8),R2S1RPPoint(-15,10,210,8)]"
CutoffRadius=32
ReversalFactor=1
Threads=16

time stack test :STCR2Z2T0S0Amodal-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${TaoDecay}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${Sigma} ${InitDist} ${CutoffRadius} ${ReversalFactor} ${Threads}  +RTS -N${Threads} -RTS"
