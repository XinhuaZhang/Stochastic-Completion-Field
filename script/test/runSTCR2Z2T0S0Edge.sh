#!/bin/bash
Size=128
Orientations=16
Scales=8
# ThetaSigma=0.3
# ScaleSigma=0.3
# MaxScale=64
# Tao=50
Trails=1000000
MaxTrails=100000
Theta0Freqs=3
ThetaFreqs=3
Scale0Freqs=3
ScaleFreqs=3
# HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat  
# HistPath=output/test/STCR2Z2T0S0Edge/${HistFileName}
NumInteration=40
WriteSourceFlag=False
EdgeFile=output/test/STCR2Z1T0Edge/avocados1.edge
Threads=16
s=4
# for s in {1,2,3}
# do 
# time stack test :STCR2Z2T0S0Edge-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${EdgeFile} ${s} ${Threads} +RTS -N${Threads} -H2g -K1g -RTS"
# # done

for MaxScale in {256,512}
do
    for Tao in {25,50,100}
    do
	for ThetaSigma in {0.2,0.3}
	do
	    for ScaleSigma in {0.2,0.3}
	    do
		HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat  
		HistPath=output/test/STCR2Z2T0S0Edge/${HistFileName}
		time stack test :STCR2Z2T0S0Edge-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${EdgeFile} ${s} ${Threads} +RTS -N${Threads} -RTS"
	    done
	done
    done
done
