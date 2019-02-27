#!/bin/bash
Size=128
Orientations=16
Sigma=0.2
Tao=100
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=100000
Threads=16

stack build Stochastic-Completion-Field:PlotGreensFunctionR2S1-exe
time stack exec -- PlotGreensFunctionR2S1-exe ${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${Threads} +RTS -N${Threads} -s -RTS
