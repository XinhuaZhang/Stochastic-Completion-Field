#!/bin/sh
Size=128
Orientations=8
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.000000001
Tao=20
Init=\(0,30,90,0\)
Trails=100000
MaxTrails=3200
SpatialFreqs=32
ThetaFreqs=10
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${SpatialFreqs}_${ThetaFreqs}_${ScaleFreqs}.dat
Threads=32

time stack test :PlotGreensFunctionSTS-test --test-arguments "${Size} ${SpatialFreqs} ${Orientations} ${ThetaSigma} ${Scales} ${ScaleSigma} ${MaxScale} ${Tao} ${Init} ${Trails} ${MaxTrails} ${ThetaFreqs} ${ScaleFreqs} ${HistFileName} ${Threads} +RTS -N${Threads} -RTS"
