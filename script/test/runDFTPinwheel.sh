#!/bin/bash
GPUID="[0]"
Size=128
NumR2Freq=128
Delta=1
DeltaRecon=1
PeriodR2=128
PeriodEnvelope=535.4916555247646
AngularFreq=5
RadialFreq=0
Sigma=-1
STD=0.5
NumBatch=8



# time stack test :DFTPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${Delta} ${DeltaRecon} ${PeriodR2} ${AngularFreq} ${RadialFreq} ${Sigma} ${PeriodEnvelope} ${STD} ${NumBatch} +RTS -N10 -RTS"
# time ./.stack-work/dist/x86_64-linux/Cabal-2.4.0.1/build/DFTPinwheel-test/DFTPinwheel-test ${GPUID} ${Size} ${NumR2Freq} ${Delta} ${DeltaRecon} ${PeriodR2} ${AngularFreq} ${RadialFreq} ${Sigma} ${PeriodEnvelope} ${STD} ${NumBatch} +RTS -N10 -RTS

theta=0
a=1
t=4
for theta in {0..71}
# for a in 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4
# for t in {0..127..4}
do
    time stack test :DFTPinwheel-test --test-arguments "${GPUID} ${Size} ${NumR2Freq} ${Delta} ${DeltaRecon} ${PeriodR2} ${AngularFreq} ${RadialFreq} ${Sigma} ${PeriodEnvelope} ${STD} ${t} ${theta} ${a} ${NumBatch} +RTS -N10 -RTS"
done

    


