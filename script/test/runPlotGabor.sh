#!/bin/bash
Freq=1
StdX=3
StdY=1
Delta=0.5
Orientation=0
Scale=1
CenterX=0
CenterY=0
NumPoints=65
Threads=2

time stack test :PlotGabor-test --test-arguments "${Freq} ${StdX} ${StdY} ${Delta} ${Orientation} ${Scale} ${CenterX} ${CenterY} ${NumPoints} +RTS -N${Threads} -RTS"

