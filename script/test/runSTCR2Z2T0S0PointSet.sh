#!/bin/bash
Size=128
Orientations=16
Scales=8
# ThetaSigma=0.3
# ScaleSigma=0.3
# MaxScale=64
# Tao=50
Trails=1000000
MaxTrails=20000
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=5
ScaleFreqs=5
# HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
# HistPath=output/test/STCR2Z2T0S0PointSet/${HistFileName}
NumInteration=20
WriteSourceFlag=False
#R=15
Threads=16

# for R in {15,20,25,30,35,40,45}
# do
# time stack test :STCR2Z2T0S0PointSet-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${R} ${Threads}  +RTS -N${Threads} -H2g -K1g -RTS"
# done

R=30

for MaxScale in {32,64,128}
do
    for Tao in {25,50,100}
    do
	for ThetaSigma in {0.0,0.1,0.2,0.3}
	do
	    for ScaleSigma in {0.0,0.1,0.2,0.3}
	    do
		HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
		HistPath=output/test/STCR2Z2T0S0PointSet/${HistFileName}
		time stack test :STCR2Z2T0S0PointSet-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${R} ${Threads}  +RTS -N${Threads} -RTS"
	    done
	done
    done
done
