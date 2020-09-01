#!/bin/bash
N=2
DeltaX=0.1
DeltaR=0.001
RadiusX=1
RadiusR=1000
Threads=2
time stack test :HankelTransform-test --test-arguments "${N} ${DeltaX} ${DeltaR} ${RadiusX} ${RadiusR} +RTS -N${Radius} -RTS"

