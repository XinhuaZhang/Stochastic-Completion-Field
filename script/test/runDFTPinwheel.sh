#!/bin/bash
GPUID="[4,5,6,7]"
Size=128
NumR2Freq=128
Delta=1
PeriodR2=128
AngularFreq=5
RadialFreq=10
Sigma=0

time stack test :DFTPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${Delta} ${PeriodR2} ${AngularFreq} ${RadialFreq} ${Sigma} +RTS -N28 -RTS"
