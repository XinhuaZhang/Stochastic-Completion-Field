#!/bin/bash
Size=128
Orientations=16
Sigma=0.1
Tao=100
InitSpeed=0
Trails=1000000
MaxTrails=100
InitOri=0
WriteFlag=True
NumInteration=15
Threshold=1
R=8
HistFileName=${Size}_${Orientations}_${Sigma}_${Tao}_${InitSpeed}_${R}.dat
HistFilePath=output/test/STCR2S1Image/${HistFileName}
ImagePath=output/test/Resize/out.png #/home/xzhang/Workspaces/Amoeba/Matlab/amoeba2D_Full/4FC/target/amoeba_1_4_00402.png
Threads=16

time stack test :STCR2S1Image-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Trails} ${MaxTrails} ${InitOri} ${InitSpeed} ${WriteFlag} ${NumInteration} ${Threshold} ${HistFilePath} ${R} ${ImagePath} ${Threads} +RTS -N${Threads} -s -RTS"
