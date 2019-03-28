#!/bin/bash
Size=256
Orientations=36
Sigma=0.1
Tao=100
Length=100
InitialScale=1
Trails=100000
MaxTrails=100000
Theta0Freqs=11
ThetaFreqs=11
# HistFileName=${Size}_${Sigma}_${Tao}_${InitialScale}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFileName=${Size}_${Sigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFilePath=output/test/STCR2Z1T0ImageSparseCode/${HistFileName}
Alpha=-3
PinwheelFlag=False
ImagePath=/home/xzhang/Workspaces/GenerativeAmoeba/output1LayerFull/Train/S1.pvp
ReconPath=/home/xzhang/Workspaces/GenerativeAmoeba/output1LayerFull/Train/S1ReconImageAmoeba.pvp
NumInteration=20
WriteSourceFlag=True
Threshold=0.05
Threads=8

time stack test :STCR2Z1T0ImageSparseCode-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${InitialScale} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${HistFilePath} ${Alpha} ${PinwheelFlag} ${ImagePath} ${ReconPath} ${NumInteration} ${WriteSourceFlag} ${Threshold} ${Threads} +RTS -N${Threads} -s -RTS"
