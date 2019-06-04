#!/bin/bash
Size=256
Orientations=16
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.000001
Tao=100
Trails=1000000 
MaxTrails=100000
Theta0Freqs=10
ThetaFreqs=10
Scale0Freqs=0
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/CannyEdgeEndPoint/${HistFileName}
NumInteration=100
WriteSourceFlag=True
CutoffRadiusEndPoint=25
CutoffRadius=64
ReversalFactor=0.01
InputImgPath=output/test/Resize/triangle_240.png
Threshold1=30
Threshold2=200
PixelDist=8
Threads=16

time stack test :CannyEdgeEndPoint-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${CutoffRadiusEndPoint} ${CutoffRadius} ${ReversalFactor} ${InputImgPath} ${Threshold1} ${Threshold2} ${PixelDist} ${Threads} +RTS -N${Threads} -RTS"
