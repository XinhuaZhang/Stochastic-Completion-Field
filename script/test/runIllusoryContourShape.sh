#!/bin/bash
GPUID="[0,1,2,3,4,5,6,7]"
Size=256
Delta=0.125 #20.03125 #0.0625
Orientations=36
Scales=36
ThetaSigma=0.0005
ScaleSigma=0.0
MaxScale=64
Cutoff=49 #49 #46 
DeltaLog=0.03125 #0.03125
Tao=9.5
Trails=96000
MaxTrails=48000
PhiFreqs=15
RhoFreqs=15
ThetaFreqs=15
ScaleFreqs=15
InitScale=0.149 #0.0078125 #0.03125
Circle="(Circle{circleNum=6,circleRadius=60})"
Shape2D="Points(0,0)2${Circle}"
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
WriteFlag=True
Iterations=20
STD=15
STDG=3
Threads=24
Suffix=64 #60
# for ThetaSigma in {0.0035,0.003,0.0025}
# do    
#     for Tao in {75,100}
#     do
	Circle="(Circle{circleNum=8,circleRadius=${Suffix}})"
	Shape2D="Points(0,0)2${Circle}"
	HistFileName="${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}.dat"
	HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
	time stack test :IllusoryContourShape-test --test-arguments "${GPUID} ${Size} ${Delta} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Cutoff} ${DeltaLog} ${Tao}  ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${InitScale} ${Shape2D} ${HistPath} ${WriteFlag} ${Iterations} ${Suffix} ${STD} ${STDG} ${Threads} +RTS -N${Threads} -RTS"
#     done
# done
