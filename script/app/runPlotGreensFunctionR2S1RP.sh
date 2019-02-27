#!/bin/bash
Size=128
Orientations=128
Scales=16
MaxScale=2
ThetaSigma=0.1
ScaleSigma=0.1
Tao=20
Step=100
Init=\(0,0,1.57,0.5,1.57,0.5\)
Trails=100000
Threads=16

stack build Stochastic-Completion-Field:PlotGreensFunctionR2S1RP-exe
time stack exec -- PlotGreensFunctionR2S1RP-exe ${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Step} ${Init} ${Trails} ${Threads} +RTS -N${Threads} -RTS
