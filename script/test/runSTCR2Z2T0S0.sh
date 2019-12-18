#!/bin/bash
Size=64
Orientations=36
Scales=32
ThetaSigma=0.01
ScaleSigma=0.0
MaxScale=32
Tao=10
Trails=100000
MaxTrails=10000
Theta0Freqs=3
ThetaFreqs=3
Scale0Freqs=3
ScaleFreqs=3
InitDist="[R2S1RPPoint(10,-10,150,1),R2S1RPPoint(-10,-10,210,1)]"
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0/${HistFileName}
Alpha=0
RadialFlag=True
Threads=32
time stack test :STCR2Z2T0S0-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${InitDist} ${HistPath} ${Alpha} ${RadialFlag} ${Threads} +RTS -N${Threads} -RTS"
Trails=0
for s in {2,4,6,8,10,12,14,16,18,20,22,24,26}
do
InitDist="[R2S1RPPoint(10,-10,150,${s}),R2S1RPPoint(-10,-10,210,${s})]"

time stack test :STCR2Z2T0S0-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${InitDist} ${HistPath} ${Alpha} ${RadialFlag} ${Threads} +RTS -N${Threads} -RTS"
done


