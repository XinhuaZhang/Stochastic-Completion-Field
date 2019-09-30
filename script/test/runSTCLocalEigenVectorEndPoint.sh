#!/bin/bash
Size=128
Orientations=16
Scales=1
ThetaSigma=0.1
ScaleSigma=0.0
MaxScale=1.000000001
Tao=100
Trails=100000
MaxTrails=100000
ThetaFreqs=10
ScaleFreqs=0
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${ThetaFreqs}_${ScaleFreqs}.dat
NumInteration=0
WriteSourceFlag=True
CutoffRadius=32
PatchNormFlag=False
PatchNormSize=160
Corner="(Corner{cornerThetaDeg0=0,cornerThetaDeg1=120,cornerLength=80})"
PacMan="(PacMan{pacManThetaDeg0=30,pacManThetaDeg1=120,pacManRadius=50})"
TJunction="(TJunction{tJunctionThetaDeg=45,tJunctionLength=50})"
Cross="(Cross{crossThetaDeg=60,crossLength=50})"
ETriangle="(ETriangle{eTriangleThetaDeg=0,eTriangleLength=100})"
IncompleteCircle="(IncompleteCircle{iCircleTheta0=0,iCircleTheta1=0,iCircleRadiau=40})"
Line="(Line{lineOrientationDeg=20,lineLength=60,lineWidth=1})"
Circle="(Circle{circleNum=8,circleRadius=35})"
Shape2D="Points(0,0)8${Circle}"
UseFFTWWisdomFlag=True
FFTWWisdomFileName="fftwWisdom.dat"
Threads=16

time stack test :STCLocalEigenVectorEndPoint-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${ThetaFreqs} ${ScaleFreqs} ${HistFileName} ${NumInteration} ${WriteSourceFlag} ${CutoffRadius} ${PatchNormFlag} ${PatchNormSize} ${Shape2D} ${UseFFTWWisdomFlag} ${FFTWWisdomFileName} ${Threads}  +RTS -N${Threads} -RTS"
