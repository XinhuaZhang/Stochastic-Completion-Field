#!/bin/bash
Size=256
Orientations=16
Sigma=0.1
Tao=100
Length=100
Init=\(-100,0,0,0,0,1.3862943611198906\)
Trails=1000000
Threads=16

stack build Stochastic-Completion-Field:PlotGreensFunctionR2S1-exe
time stack exec -- PlotGreensFunctionR2S1-exe ${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${Threads} +RTS -N${Threads} -s -RTS
