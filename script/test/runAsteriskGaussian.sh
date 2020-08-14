#!/bin/bash
R2Freq=256
ThetaFreq=20
Points=256
Theta=72
X=1
Y=0
Period=256
STD1=-1
STD2=10 #7.295125 #0.456
time stack test :AsteriskGaussian-test --test-arguments "${R2Freq} ${ThetaFreq} ${Points} ${Theta} ${X} ${Y} ${Period} ${STD1} ${STD2} +RTS -N20 -RTS"

