#!/bin/bash
GPUID="[7]"
Size=256
Delta=0.2 #0.03125 #0.0625
Orientations=36
Scales=1
ThetaSigma=0.04
ScaleSigma=0.0
MaxScale=32
Cutoff=46 #49 #46 
DeltaLog=1
Tao=40
Trails=480000
MaxTrails=48000
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=0
InitScale=0.03125 #0.0078125 #0.03125
Circle="(Circle{circleNum=6,circleRadius=60})"
Shape2D="Points(0,0)2${Circle}"
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
WriteFlag=True
Iterations=50
STD=10
Threads=24
Suffix=60 #60
# for Cutoff in {24,128,96,64,48,32}
# do    
    # for Suffix in {30,60}
    # do
	Circle="(Circle{circleNum=8,circleRadius=${Suffix}})"
	Shape2D="Points(0,0)2${Circle}"
	time stack test :IllusoryContourShape-test --test-arguments "${GPUID} ${Size} ${Delta} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Cutoff} ${DeltaLog} ${Tao}  ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitScale} ${Shape2D} ${HistPath} ${WriteFlag} ${Iterations} ${Suffix} ${STD} ${Threads} +RTS -N${Threads} -RTS"
    # done
# done
