#!/bin/bash
Size=512
Theta0Freqs=5
Alpha=-2
IDX=\(328,328\)
OrientationSamples=360
Shape="Line{lineOrientationDeg=60.0,lineLength=200,lineWidth=10}"
# Shape="PacMan{pacManThetaDeg0=0.0,pacManThetaDeg1=90.0,pacManRadius=100}"
WhiteBackgroundFlag=False
Threads=16

time stack test :OrientationAnalysisShape-test --test-arguments "${Size} ${Orientations} ${Theta0Freqs} ${Alpha} ${IDX} ${OrientationSamples} ${Shape} ${WhiteBackgroundFlag} ${Threads} +RTS -N${Threads} -RTS"
