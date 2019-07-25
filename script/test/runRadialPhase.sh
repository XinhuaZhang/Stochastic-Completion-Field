#!/bin/bash
Size=64
ThetaSigma=0.1
ScaleSigma=0.1
MaxScale=16
MaxScale1=16
Tao=100
Trails=100000
MaxTrails=100000
Theta0Freqs=3
ThetaFreqs=3
Scale0Freqs=3
ScaleFreqs=3
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/RadialPhase/${HistFileName}
Alpha=-2
DeltaTheta=0
Threads=16

time stack test :RadialPhase-test --test-arguments "${Size} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${MaxScale1} ${Tao} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${Alpha} ${DeltaTheta} ${Threads} +RTS -N${Threads} -RTS"
