#!/bin/bash
Size=128
Orientations=64
Sigma=0.1
Tao=100
Length=100
Init=\(0,0,0,0,0,1\)
Trails=100000
Theta0Freqs=7
ThetaFreqs=7
InitDist="[R2S1RPPoint(15,20,150,1),R2S1RPPoint(-15,20,30,1)]"
Threads=16

time stack test :STCR2Z1T0-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${Theta0Freqs} ${ThetaFreqs} ${InitDist} ${Threads} +RTS -N${Threads} -s -RTS"
