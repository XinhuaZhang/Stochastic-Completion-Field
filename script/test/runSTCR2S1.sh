#!/bin/bash
Size=128
Orientations=36
Sigma=0.1
Tao=100
Trails=480000   #Total number of trails
MaxTrails=48000  #The max number of trails run in one batch
InitDist="[R2S1RPPoint(25,0,150,0),R2S1RPPoint(-25,0,210,0)]" #do NOT add any space in it 
InitOri=0
InitSpeed=0
R=64
Threads=48

time stack test :STCR2S1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitDist} ${InitOri} ${InitSpeed} ${R} ${Threads} +RTS -N${Threads} -RTS"
