#!/bin/bash
Size=32
Orientations=16
ThetaSigma=0.05
ScaleSigma=0.1
MaxScale=8
Tao=100
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=100000
MaxTrails=50000
Theta0Freqs=2
ThetaFreqs=2
Scale0Freqs=2
ScaleFreqs=2
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
Alpha=-1
Threads=16

stack build Stochastic-Completion-Field:PlotGreensFunctionR2Z2T0S0-exe
time stack exec -- PlotGreensFunctionR2Z2T0S0-exe ${Size} ${Orientations} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistFileName} ${Alpha} ${Threads} +RTS -N${Threads} -RTS
