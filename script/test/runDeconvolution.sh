#!/bin/bash
Size=64
Orientations=72
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.0000001
Tao=100
Trails=100000
MaxTrails=100000
Theta0Freqs=10
ThetaFreqs=10
Scale0Freqs=0
ScaleFreqs=0
InitDist="[R2S1RPPoint(8,8,120,1),R2S1RPPoint(0,0,315,1)]" #R2S1RPPoint(-15,10,210,8),R2S1RPPoint(0,0,315,1)
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/Deconvolution/${HistFileName}
Alpha=-2
PinwheelFlag=True
Threads=16
# s=1
time stack test :Deconvolution-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${InitDist} ${HistPath} ${Alpha} ${PinwheelFlag} ${Threads} +RTS -N${Threads} -H2g -K1g -RTS"

