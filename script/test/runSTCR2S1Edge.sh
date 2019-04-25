#!/bin/bash
Size=512
Orientations=16
Sigma=0.1
Tao=100
InitSpeed=1
Trails=1000000
MaxTrails=100
InitOri=0
WriteFlag=True
NumInteration=5
Threshold=1
HistFileName=${Size}_${Orientations}_${Sigma}_${Tao}_${InitSpeed}.dat
HistFilePath=output/test/STCR2S1Edge/${HistFileName}
EdgeFile=output/test/STCR2Z1T0Edge/avocados1.edge
Threads=16

time stack test :STCR2S1Edge-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitOri} ${InitSpeed} ${WriteFlag} ${NumInteration} ${Threshold} ${HistFilePath} ${EdgeFile} ${Threads} +RTS -N${Threads} -s -RTS"
