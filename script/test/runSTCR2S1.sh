#!/bin/bash
Size=128
Orientations=72
Sigma=0.0005
Tao=9.5
Trails=4800000   #Total number of trails
MaxTrails=480000  #The max number of trails run in one batch
InitDist="[R2S1RPPoint(-60,0,0,0),R2S1RPPoint(30,0,135,0)]" #do NOT add any space in it 
InitOri=0
InitSpeed=0.149
R=61
Delta=1
Location="(25,4)"
Threads=24

time stack test :STCR2S1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitDist} ${InitOri} ${InitSpeed} ${R} ${Location} ${Delta} ${Threads} +RTS -N${Threads} -RTS"
