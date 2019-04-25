#!/bin/bash
Size=64
Orientations=36
Scales=32
ThetaSigma=0.1
ScaleSigma=0.2
MaxScale=16
Tao=20
Trails=100000
MaxTrails=2000
Theta0Freqs=5
ThetaFreqs=5
Scale0Freqs=5
ScaleFreqs=5
# InitDist="[R2S1RPPoint(15,10,150,8),R2S1RPPoint(-15,10,210,8)]"
HistFileName=${Size}_${ThetaSigma}_${ScaleSigma}_${MaxScale}_${Tao}_${Theta0Freqs}_${ThetaFreqs}_${Scale0Freqs}_${ScaleFreqs}.dat
HistPath=output/test/STCR2Z2T0S0/${HistFileName}
Alpha=-2
PinwheelFlag=True
Threads=16
# s=1
for s in {1,2,4,6,8,10,12,14,16,18,20,22,24,26}
do
InitDist="[R2S1RPPoint(15,10,150,${s}),R2S1RPPoint(-15,10,210,${s})]"

time stack test :STCR2Z2T0S0-test --test-arguments "${Size} ${Orientations} ${Scales} ${ThetaSigma} ${ScaleSigma} ${MaxScale} ${Tao}  ${Trails} ${MaxTrails} ${Theta0Freqs} ${ThetaFreqs} ${Scale0Freqs} ${ScaleFreqs} ${InitDist} ${HistPath} ${Alpha} ${PinwheelFlag} ${Threads} +RTS -N${Threads} -H2g -K1g -RTS"
done
