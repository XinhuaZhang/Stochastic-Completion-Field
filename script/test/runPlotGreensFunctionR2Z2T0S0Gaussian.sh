#!/bin/bash
Size=64
Orientations=16
ThetaSigma=0.08
ScaleSigma=0.0
MaxScale=1.00000001
Tao=15
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=1000000
MaxTrails=100000
Theta0Freqs=3
ThetaFreqs=3
Scale0Freqs=0
ScaleFreqs=0
Sigma=4
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}_${Sigma}.dat
Alpha=-2
Threads=16

time stack test :PlotGreensFunctionR2Z2T0S0Gaussian-test --test-arguments "${Size} ${Orientations} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistFileName} ${Alpha} ${Sigma} ${Threads} +RTS -N${Threads} -RTS"
