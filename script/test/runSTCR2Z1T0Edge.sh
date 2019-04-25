#!/bin/bash
Size=512
Orientations=36
Sigma=0.1
Tao=100
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
NumInteration=50
WriteSourceFlag=True
EdgeFile=output/test/STCR2Z1T0Edge/avocados1.edge
Threads=16

time stack test :STCR2Z1T0Edge-test --test-arguments "${Size} ${Orientations} ${Sigma} ${Tao} ${Length} ${InitialScale} ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${HistFilePath} ${PinwheelFlag} ${NumInteration} ${WriteSourceFlag} ${EdgeFile} ${Threads} +RTS -N${Threads} -s -RTS"
