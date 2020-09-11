#!/bin/bash
R2Freq=256
ThetaFreq=41
RFreq=41
Points=256
Theta=72
R=150
X=0.5 #1.831563888873418e-2 #0.1353352832366127 #7.38905609893065
Y=0
Period=256
STD1=-1
STD2=0.5
#7.295125 #0.456
STDR2=0.25
time stack test :AsteriskGaussian-test --test-arguments "${R2Freq} ${ThetaFreq} ${RFreq} ${Points} ${Theta} ${R} ${X} ${Y} ${Period} ${STD1} ${STD2} ${STDR2} +RTS -N28 -RTS"

