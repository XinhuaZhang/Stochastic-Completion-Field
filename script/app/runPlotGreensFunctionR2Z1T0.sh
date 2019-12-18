#!/bin/bash
Size=128
Orientations=16
Sigma=0.1
Tao=100
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=100000
MaxTrails=100000
Theta0Freqs=[-2,-1,0,1,2]
ThetaFreqs=[-2,-1,0,1,2]
Threads=16

stack build Stochastic-Completion-Field:PlotGreensFunctionR2Z1T0-exe
time stack exec -- PlotGreensFunctionR2Z1T0-exe ${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Threads} +RTS -N${Threads} -RTS
