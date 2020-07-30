#!/bin/bash
R2Freq=256
ThetaFreq=101
Points=256
Theta=72
X=0.5
Y=0.25
Period=256
STD1=-1
STD2=7.295125 #0.456
time stack test :AsteriskGaussian-test --test-arguments "${R2Freq} ${ThetaFreq} ${Points} ${Theta} ${X} ${Y} ${Period} ${STD1} ${STD2} +RTS -N28 -RTS"

