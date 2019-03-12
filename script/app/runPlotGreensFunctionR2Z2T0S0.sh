#!/bin/bash
Size=128
Orientations=16
ThetaSigma=0.1
ScaleSigma=0.05
MaxScale=91
Tao=100
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=100000
MaxTrails=480
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=5
ScaleFreqs=5
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
Alpha=-1
Threads=16

stack build Stochastic-Completion-Field:PlotGreensFunctionR2Z2T0S0-exe
time stack exec -- PlotGreensFunctionR2Z2T0S0-exe ${Size} ${Orientations} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistFileName} ${Alpha} ${Threads} +RTS -N${Threads} -s -RTS
