#!/bin/bash
Size=128
Orientations=36
Sigma=0.1
Tao=20
Length=100
InitialScale=1
Trails=1000000
MaxTrails=100000
Theta0Freqs=11
ThetaFreqs=11
# HistFileName=${Size}_${Sigma}_${Tao}_${InitialScale}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFileName=${Size}_${Sigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}.dat
HistFilePath=output/test/STCR2Z1T0PointSet/${HistFileName}
PinwheelFlag=False
NumInteration=40
WriteSourceFlag=True
Threads=8

time stack test :STCR2Z1T0PointSet-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${InitialScale} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${HistFilePath} ${PinwheelFlag} ${NumInteration} ${WriteSourceFlag} ${Threads} +RTS -N${Threads} -s -RTS"
