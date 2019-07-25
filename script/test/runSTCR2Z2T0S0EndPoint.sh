#!/bin/bash
Size=128
Orientations=16
Scales=16
ThetaSigma=0.1
ScaleSigma=0.1
MaxScale=32 #1.000000001
Tao=500
Trails=1000000
MaxTrails=10000
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=5
ScaleFreqs=5
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0EndPoint/${HistFileName}
NumInteration=10
WriteSourceFlag=True
CutoffRadiusEndPoint=32
CutoffRadius=32
ReversalFactor=0.0
PatchNormFlag=False
PatchNormSize=160
ApproximatedEigenValue=0.0
Corner="(Corner{cornerThetaDeg0=0,cornerThetaDeg1=90,cornerLength=80})"
PacMan="(PacMan{pacManThetaDeg0=30,pacManThetaDeg1=120,pacManRadius=50})"
TJunction="(TJunction{tJunctionThetaDeg=45,tJunctionLength=50})"
Cross="(Cross{crossThetaDeg=60,crossLength=50})"
ETriangle="(ETriangle{eTriangleThetaDeg=0,eTriangleLength=100})"
IncompleteCircle="(IncompleteCircle{iCircleTheta0=0,iCircleTheta1=0,iCircleRadiau=31})"
Line="(Line{lineOrientationDeg=20,lineLength=90,lineWidth=1})"
Shape2D="Points(0,0)21${IncompleteCircle}"

Threads=16
# for R in {30,35}
# do
time stack test :STCR2Z2T0S0EndPoint-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${CutoffRadiusEndPoint} ${CutoffRadius} ${ReversalFactor} 60 ${PatchNormFlag} ${PatchNormSize} ${ApproximatedEigenValue} ${Shape2D} ${Threads}  +RTS -N${Threads} -H2g -K1g -RTS"
# done

# R=30

# for MaxScale in {6,10,14,18,22,26,30}
# do
    # for Tao in {50,100}
    # do
	# for ThetaSigma in {0.1,0.125,0.075}
	# do
	    # for ScaleSigma in {0.25,0.3}
	    # do
	    # 	HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
	    # 	HistPath=output/test/STCR2Z2T0S0EndPoint/${HistFileName}
            #     # rm ${HistPath}
	    # 	time stack test :STCR2Z2T0S0EndPoint-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${R} ${Threads}  +RTS -N${Threads} -RTS"
	    # done
	# done
    # done
# done
