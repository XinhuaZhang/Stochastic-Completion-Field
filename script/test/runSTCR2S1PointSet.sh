#!/bin/bash
Size=256
Orientations=72
Sigma=0.0025 #0.0005 #0.00447213595499958 #0.0005
Tao=75 # 9.5
InitSpeed=1 #0.149
Trails=1000000
MaxTrails=100000
InitOri=0
WriteFlag=True
NumInteration=5
Threshold=1
HistFileName=${Size}_${Orientations}_${Sigma}_${Tao}_${InitSpeed}.dat
HistFilePath=output/test/STCR2S1PointSet/${HistFileName}
R=64
Delta=1 #0.125
STD=0.125
Threads=24

time stack test :STCR2S1PointSet-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitOri} ${InitSpeed} ${WriteFlag} ${NumInteration} ${Threshold} ${HistFilePath} ${R} ${Delta} ${STD} ${Threads} +RTS -N${Threads} -RTS"
