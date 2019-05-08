#!/bin/bash
Size=64
Orientations=16
Sigma=0.1
Tao=100
Trails=1000000
MaxTrails=10000
InitDist="[R2S1RPPoint(0,0,0,0),R2S1RPPoint(-20,10,160,0)]"
InitOri=0
InitSpeed=0
R=6
Threads=16

time stack test :STCR2S1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitDist} ${InitOri} ${InitSpeed} ${R} ${Threads} +RTS -N${Threads} -s -RTS"
