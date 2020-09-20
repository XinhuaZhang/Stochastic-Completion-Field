#!/bin/bash
R2Freq=256
ThetaFreq=41
RFreq=41
Points=256
Theta=72
R=150
AngularFreq=0
RadialFreq=0
STD=1.4
Sigma=1
# time stack test :GaussianPinwheel-test --test-arguments "${R2Freq} ${ThetaFreq} ${RFreq} ${Points} ${Theta} ${R} ${AngularFreq} ${RadialFreq} ${STD} ${Sigma} +RTS -N14 -RTS"

time .stack-work/dist/x86_64-linux/Cabal-2.4.0.1/build/GaussianPinwheel-test/GaussianPinwheel-test ${R2Freq} ${ThetaFreq} ${RFreq} ${Points} ${Theta} ${R} ${AngularFreq} ${RadialFreq} ${STD} ${Sigma} +RTS -N14 -RTS
