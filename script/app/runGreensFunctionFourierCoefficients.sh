#!/bin/bash
GPU="[7]" #"[0,1,2,3,4,5,6,7]"
ThetaSigma=0.025
ScaleSigma=0.0
MaxScale=64
DeltaLog=0.1
Tao=100
Trails=240000 #768000 #96000
MaxTrails=48000 #38400 #4800
PhiFreqs=15
RhoFreqs=10
ThetaFreqs=15
ScaleFreqs=10
InitScale=0.1
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}_${InitScale}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
Threads=48
stack build Stochastic-Completion-Field:GreensFunctionFourierCoefficients-exe
time stack exec -- GreensFunctionFourierCoefficients-exe ${GPU} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${DeltaLog} ${Tao} ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs}  ${HistPath} ${InitScale} ${Threads} +RTS -N${Threads} -RTS
