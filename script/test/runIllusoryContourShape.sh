#!/bin/bash
Size=256
Delta=1
Orientations=8
Scales=8
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=128
DeltaLog=0.1
Tao=100
Trails=14400
MaxTrails=4800
PhiFreqs=15
RhoFreqs=10
ThetaFreqs=15
ScaleFreqs=0
Circle="(Circle{circleNum=8,circleRadius=30})"
Shape2D="Points(0,0)2${Circle}"
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
WriteFlag=True
Iterations=10
Threads=48
time stack test :IllusoryContourShape-test --test-arguments "${Size} ${Delta} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${DeltaLog} ${Tao}  ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${Shape2D} ${HistPath} ${WriteFlag} ${Iterations} ${Threads} +RTS -N${Threads} -RTS"
