#!/bin/bash
GPUID="[2,3,4,5,6,7]"
Size=256 #2047
SizeRecon=256
Delta=1 #0.03125
DeltaRecon=1
NumR2Freq=256
PeriodR2=256
PhiFreqs=0
RhoFreqs=0
STDR2=-0.5 #0.03125
STD=15
S=1.5
NumBatchR2=2
NumBatchR2Freqs=1
Threads=28

time stack test :EnvelopeTest-test --test-arguments "${GPUID} ${Size} ${Delta} ${SizeRecon} ${DeltaRecon}  ${NumR2Freq} ${PeriodR2} ${PhiFreqs} ${RhoFreqs} ${STDR2} ${STD} ${NumBatchR2} ${NumBatchR2Freqs} ${S} ${Threads} +RTS -N${Threads} -RTS"

