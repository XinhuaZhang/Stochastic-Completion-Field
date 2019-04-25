#!/bin/bash
Size=64
Orientations=64
Sigma=0.1
Tao=20
Length=100
Init=\(0,0,0,0,0,1\)
Trails=100000
MaxTrails=100000
Theta0Freqs=10
ThetaFreqs=10
InitDist="[R2S1RPPoint(-15,10,30,1),R2S1RPPoint(15,10,330,1)]"
HistFileName=${Size}_${Sigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFilePath=output/test/STCR2Z1T0/${HistFileName}
Alpha=-2
PinwheelFlag=False
Threads=16

time stack test :STCR2Z1T0-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${Init} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${InitDist} ${HistFilePath} ${Alpha} ${PinwheelFlag} ${Threads} +RTS -N${Threads} -s -RTS"
