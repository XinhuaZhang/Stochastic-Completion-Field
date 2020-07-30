#!/bin/bash
GPUID="[7]"
Size=128
NumR2Freq=256
PeriodR2=128
PhiFreqs=10
RhoFreqs=10
ThetaFreqs=10
ScaleFreqs=0
STD=0.0078125
Radius=0.015625
S=1
Idx="(5,0)"
Batch=16
A=1
B=50
Delta=1
RadialFreq=0
Threads=28

time stack test :HollowPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${STD} ${Radius} ${S} ${Idx} ${Batch} ${A} ${B} ${Delta} ${RadialFreq} ${Threads} +RTS -N${Threads} -RTS"

