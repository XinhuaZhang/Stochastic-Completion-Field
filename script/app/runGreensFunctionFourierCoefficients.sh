#!/bin/bash
ThetaSigma=0.1
ScaleSigma=0.05
MaxScale=64
DeltaLog=0.1
Tao=10
Trails=24000
MaxTrails=4800
PhiFreqs=5
RhoFreqs=3
ThetaFreqs=5
ScaleFreqs=3
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
Threads=48
stack build Stochastic-Completion-Field:GreensFunctionFourierCoefficients-exe
time stack exec -- GreensFunctionFourierCoefficients-exe ${Size} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${DeltaLog} ${Tao} ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs}  ${HistPath} ${Threads} +RTS -N${Threads} -RTS
