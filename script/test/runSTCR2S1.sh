#!/bin/bash
Size=64
Orientations=16
Sigma=0.1
Tao=100
Trails=1000000   #Total number of trails
MaxTrails=10000  #The max number of trails run in one batch
InitDist="[R2S1RPPoint(-10,10,45,0),R2S1RPPoint(10,10,315,0)]" #do NOT add any space in it 
InitOri=0
InitSpeed=0
R=32
Threads=16

time stack test :STCR2S1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitDist} ${InitOri} ${InitSpeed} ${R} ${Threads} +RTS -N${Threads} -RTS"
