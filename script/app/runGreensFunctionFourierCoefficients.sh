#!/bin/bash
GPU=7
ThetaSigma=0.1
ScaleSigma=0.1
MaxScale=128
DeltaLog=0.1
Tao=20
Trails=96000
MaxTrails=4800
PhiFreqs=10
RhoFreqs=5
ThetaFreqs=10
ScaleFreqs=5
HistFileName=${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${PhiFreqs}_${RhoFreqs}_${ThetaFreqs}_${ScaleFreqs}.dat
HistPath=output/app/GreensFunctionFourierCoefficients/${HistFileName}
Threads=48
stack build Stochastic-Completion-Field:GreensFunctionFourierCoefficients-exe
time stack exec -- GreensFunctionFourierCoefficients-exe ${GPU} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${DeltaLog} ${Tao} ${Trails} ${MaxTrails} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs}  ${HistPath} ${Threads} +RTS -N${Threads} -RTS
