#!/bin/bash
Size=128
Orientations=16
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.000000001
Tao=100
Trails=100000
MaxTrails=100000
ThetaFreqs=20
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${ThetaFreqs}_${ScaleFreqs}.dat
Threads=16

time stack test :PlotLocalEigenVector-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Trails} ${MaxTrails} ${ThetaFreqs} ${ScaleFreqs} ${HistFileName} ${Threads}  +RTS -N${Threads}  -RTS"
