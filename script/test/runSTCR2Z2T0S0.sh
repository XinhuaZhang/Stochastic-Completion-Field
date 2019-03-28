#!/bin/bash
Size=64
Orientations=64
Scales=32
ThetaSigma=0.1
ScaleSigma=0.05
MaxScale=45
Tao=100
Length=100
Init=\(0,0,0,0,0,1\)
Trails=100000
MaxTrails=1000
Theta0Freqs=11
ThetaFreqs=11
Scale0Freqs=11
ScaleFreqs=11
InitDist="[R2S1RPPoint(15,10,30,1),R2S1RPPoint(-15,10,330,1)]"
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0/${HistFileName}
Alpha=-2
PinwheelFlag=True
Threads=8

time stack test :STCR2Z2T0S0-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${InitDist} ${HistPath} ${Alpha} ${PinwheelFlag} ${Threads} +RTS -N${Threads} -s -RTS"
