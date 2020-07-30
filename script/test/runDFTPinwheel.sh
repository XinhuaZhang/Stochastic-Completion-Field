#!/bin/bash
GPUID="[4,5,6,7]"
Size=512
NumR2Freq=1024
Delta=1
PeriodR2=512
PeriodEnvelope=256
AngularFreq=-5
RadialFreq=0
Sigma=1
NumBatch=128

time stack test :DFTPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${Delta} ${PeriodR2} ${AngularFreq} ${RadialFreq} ${Sigma} ${PeriodEnvelope} ${NumBatch} +RTS -N28 -RTS"
