#!/bin/bash
Size=128
Orientations=16
Scales=16
ThetaSigma=0.075
ScaleSigma=0.25
MaxScale=32
TaoDecay=100
TaoReversal=200
TaoCorner=50
Trails=1000000
MaxTrails=100000
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=3
ScaleFreqs=3
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${TaoDecay}_${TaoReversal}_${TaoCorner}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0KoffkaCross/${HistFileName}
NumInteration=20
WriteSourceFlag=False
D=20
W=12
Sigma=0.1
Theta=0
Threads=16

for W in {4,6,8,10,12}
do
time stack test :STCR2Z2T0S0KoffkaCross-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${TaoDecay} ${TaoReversal} ${TaoCorner}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${D} ${W} ${Sigma} ${Theta} ${Threads}  +RTS -N${Threads} -RTS"
done

# R=30

# for MaxScale in {32,64,128}
# do
#     for Tao in {25,50,100}
#     do
# 	for ThetaSigma in {0.0,0.1,0.2,0.3}
# 	do
# 	    for ScaleSigma in {0.0,0.1,0.2,0.3}
# 	    do
# 		HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
# 		HistPath=output/test/STCR2Z2T0S0KoffkaCross/${HistFileName}
# 		time stack test :STCR2Z2T0S0KoffkaCross-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${R} ${Threads}  +RTS -N${Threads} -RTS"
# 	    done
# 	done
#     done
# done
