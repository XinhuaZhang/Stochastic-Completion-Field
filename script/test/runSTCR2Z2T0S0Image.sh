#!/bin/bash
Size=128
Orientations=16
Scales=16
ThetaSigma=0.1
ScaleSigma=0.3
MaxScale=8
Tao=100
Trails=1000000
MaxTrails=100000
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=3
ScaleFreqs=3
# HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
# HistPath=output/test/STCR2Z2T0S0Image/${HistFileName}
NumInteration=20
WriteSourceFlag=False
ImagePath=output/test/Resize/out.png #/home/xzhang/Workspaces/Amoeba/Matlab/amoeba2D_Full/4FC/target/amoeba_1_4_00402.png
CutoffR=4
ReversalFactor=0.001
Threads=16

# for MaxScale in {32,64,128}
# do
#     for Tao in {25,50,100}
#     do
# 	for ThetaSigma in {0.2,0.3}
# 	do
# 	    for ScaleSigma in {0.1,0.2,0.3}

# 	    do
		HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
		HistPath=output/test/STCR2Z2T0S0Image/${HistFileName}
		time stack test :STCR2Z2T0S0Image-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${HistPath} ${NumInteration} ${WriteSourceFlag} ${ImagePath} ${CutoffR} ${ReversalFactor} ${Threads}  +RTS -N${Threads} -H2g -K1g -RTS"
# 	    done
# 	done
#     done
# done

