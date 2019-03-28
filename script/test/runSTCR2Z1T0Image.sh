#!/bin/bash
Size=256
Orientations=36
Sigma=0.1
Tao=20
Length=100
InitialScale=1
Trails=100000
MaxTrails=100000
Theta0Freqs=5
ThetaFreqs=5
# HistFileName=${Size}_${Sigma}_${Tao}_${InitialScale}_${Theta0Freqs}_${ThetaFreqs}.dat
# HistFileName=${Size}_${Sigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFileName=test.dat
HistFilePath=output/test/STCR2Z1T0Image/${HistFileName}
Alpha=-3
PinwheelFlag=False
ImagePath=/home/xzhang/Downloads/196073.jpg #/home/xzhang/Workspaces/Amoeba/Matlab/amoeba2D_Full/4FC/target/amoeba_1_4_00402.png #output/test/Resize/out.png #output/test/STCR2Z1T0Image/Completion_normalized.png
NumInteration=40
WriteSourceFlag=True
Threads=16

time stack test :STCR2Z1T0Image-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${InitialScale} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${HistFilePath} ${Alpha} ${PinwheelFlag} ${ImagePath} ${NumInteration} ${WriteSourceFlag} ${Threads} +RTS -N${Threads} -s -RTS"
