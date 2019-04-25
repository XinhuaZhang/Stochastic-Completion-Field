#!/bin/bash
Size=128
Orientations=36
Sigma=0.1
Tao=100
Trails=100000
MaxTrails=10000
InitDist="[R2S1RPPoint(15,10,210,1),R2S1RPPoint(-15,10,150,1)]"
InitOri=0
InitSpeed=8
Threads=16

time stack test :STCR2S1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitDist} ${InitOri} ${InitSpeed} ${Threads} +RTS -N${Threads} -s -RTS"
