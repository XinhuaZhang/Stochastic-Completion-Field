#!/bin/bash
Size=512
Orientations=16
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.00000001
Tao=100
Trails=1000000
MaxTrails=100000
ThetaFreqs=10
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${ThetaFreqs}_${ScaleFreqs}.dat  
HistPath=output/test/STCR2Z2T0S0Edge/${HistFileName}
NumInteration=40
WriteSourceFlag=True
EdgeFile=output/test/STCR2Z2T0S0Edge/avocados1.edge
NumNoisePoint=100
S=1
UseFFTWWisdomFlag=True
FFTWWisdomFileName="fftwWisdom.dat"
Threads=32

# for s in {1,2,3}
# do 
time stack test :STCR2Z2T0S0Edge-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${ThetaFreqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${EdgeFile} ${NumNoisePoint} ${S} ${UseFFTWWisdomFlag} ${FFTWWisdomFileName} ${Threads} +RTS -N${Threads} -RTS" --ghc-options -j8
# # done

# for MaxScale in {256,512}
# do
#     for Tao in {25,50,100}
#     do
# 	for ThetaSigma in {0.2,0.3}
# 	do
# 	    for ScaleSigma in {0.2,0.3}
# 	    do
# 		HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat  
# 		HistPath=output/test/STCR2Z2T0S0Edge/${HistFileName}
# 		time stack test :STCR2Z2T0S0Edge-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${EdgeFile} ${s} ${Threads} +RTS -N${Threads} -RTS"
# 	    done
# 	done
#     done
# done
