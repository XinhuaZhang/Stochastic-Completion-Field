#!/bin/bash
Size=128
Orientations=64
Sigma=0.1
Tao=100
Length=100
Init=\(0,0,0,0,0,1\)
Trails=100000
MaxTrails=100000
Theta0Freqs=5
ThetaFreqs=5
InitDist="[R2S1RPPoint(15,20,180,1),R2S1RPPoint(-15,20,180,1)]"
HistFileName=${Size}_${Sigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFilePath=output/test/STCR2Z1T0/${HistFileName}
Alpha=-2
PinwheelFlag=False
Threads=16

time stack test :STCR2Z1T0_1-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${InitDist} ${HistFilePath} ${Alpha} ${PinwheelFlag} ${Threads} +RTS -N${Threads} -s -RTS"
