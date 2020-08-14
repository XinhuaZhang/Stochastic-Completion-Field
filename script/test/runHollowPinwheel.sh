#!/bin/bash
GPUID="[0]"
Size=256
NumR2Freq=256
PeriodR2=256
PhiFreqs=20
RhoFreqs=20
ThetaFreqs=20
ScaleFreqs=0
STD=0.0078125
Radius=8
S=1
Idx="(0,5)"
Batch=16
A=1
B=50
Delta=1
RadialFreq=0
Threads=20

time stack test :HollowPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${ThetaFreqs} ${ScaleFreqs} ${STD} ${Radius} ${S} ${Idx} ${Batch} ${A} ${B} ${Delta} ${RadialFreq} ${Threads} +RTS -N${Threads} -RTS"

