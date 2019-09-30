#!/bin/bash
Size=256
Orientations=16
Scales=16
ThetaSigma=0.1
ScaleSigma=0.1
MaxScale=48 #1.0000001
Tao=100
Trails=1000000 
MaxTrails=100000
Theta0Freqs=8
ThetaFreqs=8
Scale0Freqs=0
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistFileNameEndPoint=${Size}_${ThetaSigma}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_EndPoint.dat
NumInteration=20
NumInterationEndPoint=100
WriteSourceFlag=True
CutoffRadiusEndPoint=16
CutoffRadius=48
ReversalFactor=0.0
InputImgPath=output/test/Resize/Ehrenstein_only.png
Threshold1=30
Threshold2=200
WriteSegmentsFlag=True
WriteEndPointFlag=True
MinimumSegmentLength=5
UseFFTWWisdomFlag=False
FFTWWisdomFileName="fftwWisdom.dat"
Threads=16

time stack test :CannyEdgeIllusoryContour-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistFileName} ${HistFileNameEndPoint} ${NumInteration} ${NumInterationEndPoint} ${WriteSourceFlag} ${CutoffRadiusEndPoint} ${CutoffRadius} ${ReversalFactor} ${InputImgPath} ${Threshold1} ${Threshold2} ${WriteSegmentsFlag} ${WriteEndPointFlag} ${MinimumSegmentLength} ${UseFFTWWisdomFlag} ${FFTWWisdomFileName} ${Threads} +RTS -N${Threads} -RTS"
