#!/bin/bash
Size=64
Orientations=8
Scales=8
ThetaSigma=0.05
ScaleSigma=0.25
MaxScale=8
Tao=100
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=1000000
MaxTrails=2000
Theta0Freqs=3
ThetaFreqs=3
Scale0Freqs=3
ScaleFreqs=3
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/R2Z2T0S0ToR2S1RPT0S0/${HistFileName}
Alpha=-1
PinwheelFlag=True
Threads=16

time stack test :R2Z2T0S0ToR2S1RPT0S0-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${Alpha} ${PinwheelFlag} ${Threads} +RTS -N${Threads} -s -RTS"
