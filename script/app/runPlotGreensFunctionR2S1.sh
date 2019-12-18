#!/bin/bash
Size=256
Orientations=16
Sigma=0.1
Tao=10
Length=100
Init=\(-100,0,0,0,0,32\)
Trails=100000
Threads=32

stack build Stochastic-Completion-Field:PlotGreensFunctionR2S1-exe
for s in {16,32}
do
    Init=\(-100,0,0,0,0,${s}\)
    time stack exec -- PlotGreensFunctionR2S1-exe ${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${Threads} +RTS -N${Threads} -s -RTS
done
