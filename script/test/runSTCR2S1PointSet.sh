#!/bin/bash
Size=64
Orientations=18
Sigma=0.1
Tao=100
InitSpeed=0
Trails=1000000
MaxTrails=100
InitOri=0
WriteFlag=True
NumInteration=50
Threshold=1
HistFileName=${Size}_${Orientations}_${Sigma}_${Tao}_${InitSpeed}.dat
HistFilePath=output/test/STCR2S1PointSet/${HistFileName}
R=6
Threads=16

time stack test :STCR2S1PointSet-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitOri} ${InitSpeed} ${WriteFlag} ${NumInteration} ${Threshold} ${HistFilePath} ${R} ${Threads} +RTS -N${Threads} -s -RTS"
