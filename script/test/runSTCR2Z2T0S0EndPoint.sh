#!/bin/bash
Size=64
Orientations=16
Scales=8
ThetaSigma=0.1
ScaleSigma=0.2
MaxScale=16 #1.000001
Tao=100
Trails=1000000 
MaxTrails=100000
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=5
ScaleFreqs=5
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0EndPoint/${HistFileName}
NumInteration=100
WriteSourceFlag=True
CutoffRadiusEndPoint=12
CutoffRadius=48
ReversalFactor=0.001
Threads=16

# for R in {30,35}
# do
time stack test :STCR2Z2T0S0EndPoint-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${CutoffRadiusEndPoint} ${CutoffRadius} ${ReversalFactor} 60 ${Threads}  +RTS -N${Threads} -H2g -K1g -RTS"
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
