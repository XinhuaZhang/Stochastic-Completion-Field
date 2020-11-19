#!/bin/bash
Size=128
Orientations=72
Sigma=0.015
Tao=10
Threshold=0.0000001
Trails=4800000   #Total number of trails
MaxTrails=480000  #The max number of trails run in one batch
InitDist="[R2S1RPPoint(-8,0,0,0),R2S1RPPoint(0,-8,90,0)]" #do NOT add any space in it 
InitOri=0
InitSpeed=1
R=61
Delta=0.125
Location="(0,0)"
Poisson=0.00001
Threads=10

time stack test :STCR2S1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Threshold} ${Trails} ${MaxTrails} ${InitDist} ${InitOri} ${InitSpeed} ${R} ${Location} ${Delta} ${Poisson} ${Threads} +RTS -N${Threads} -RTS"
