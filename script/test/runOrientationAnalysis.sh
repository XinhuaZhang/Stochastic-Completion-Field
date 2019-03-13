#!/bin/bash
Size=256
Orientations=100
Sigma=0.1
Tao=20
Length=100
Init=\(0,0,0,0,0,1\)
Trails=100000
MaxTrails=100000
Theta0Freqs=3
ThetaFreqs=3
HistFileName=${Size}_${Sigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFilePath=output/test/STCR2Z1T0Image/${HistFileName}
Alpha=-3
PinwheelFlag=False
ImagePath=/home/xzhang/Workspaces/Amoeba/Matlab/amoeba2D_Full/4FC/target/amoeba_1_4_00201.png #output/test/Resize/out.png
IDX=\(64,64\)
OrientationSamples=100
Threads=16

time stack test :OrientationAnalysis-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${HistFilePath} ${Alpha} ${PinwheelFlag} ${ImagePath} ${IDX} ${OrientationSamples} ${Threads} +RTS -N${Threads} -s -RTS"
