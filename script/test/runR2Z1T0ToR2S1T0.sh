#!/bin/bash
Size=128
Orientations=16
Sigma=0.1
Tao=100
Length=100
Init=\(0,0,0,0,1.57,1\)
Trails=10000
Theta0Freqs=[-3,-2,-1,0,1,2,3]
ThetaFreqs=[-3,-2,-1,0,1,2,3]
Threads=16

time stack test :R2Z1T0ToR2S1T0-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${Theta0Freqs} ${ThetaFreqs} ${Threads} +RTS -N${Threads} -s -RTS"
