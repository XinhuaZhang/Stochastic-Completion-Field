#!/bin/bash
Size=128
Orientations=16
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.000000001
Tao=100
Trails=1000000
MaxTrails=100000
ThetaFreqs=10
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${ThetaFreqs}_${ScaleFreqs}.dat
CutoffRadius=64
InitDist="[R2S1RPPoint(15,10,150,1),R2S1RPPoint(-15,10,30,1)]"
UseFFTWWisdomFlag=True
FFTWWisdomFileName="fftwWisdom.dat"
Threads=16

time stack test :STCLocalEigenVector-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao} ${Trails} ${MaxTrails} ${ThetaFreqs} ${ScaleFreqs} ${HistFileName} ${CutoffRadius} ${InitDist} ${UseFFTWWisdomFlag} ${FFTWWisdomFileName} ${Threads} +RTS -N${Threads} -RTS"
