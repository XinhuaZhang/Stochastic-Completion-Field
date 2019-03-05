#!/bin/bash
Size=64
Orientations=16
Scales=16
ThetaSigma=0.1
ScaleSigma=0.05
MaxScale=45
Tao=100
Length=100
Init=\(0,0,0,0,0,1\)
Trails=100000
MaxTrails=50000
Theta0Freqs=3
ThetaFreqs=3
Scale0Freqs=3
ScaleFreqs=3
InitDist="[R2S1RPPoint(15,20,150,1),R2S1RPPoint(-15,20,30,1)]"
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/app/PlotGreensFunctionR2Z2T0S0/${HistFileName}
Threads=16

time stack test :STCR2Z2T0S0-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${InitDist} ${HistPath} ${Threads} +RTS -N${Threads} -s -RTS"
