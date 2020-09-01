#!/bin/bash
GPUID="[0]"
Size=1024
NumR2Freq=1024
Delta=1
DeltaRecon=1
PeriodR2=1024
PeriodEnvelope=1024
AngularFreq=0
RadialFreq=0
Sigma=1
STD=0.5
NumBatch=128

time stack test :DFTPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${Delta} ${DeltaRecon} ${PeriodR2} ${AngularFreq} ${RadialFreq} ${Sigma} ${PeriodEnvelope} ${STD} ${NumBatch} +RTS -N20 -RTS"
