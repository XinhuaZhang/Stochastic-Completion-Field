#!/bin/bash
InputPath=/home/xzhang/Workspaces/ComputerVision/Application/ECCV2018/Presentation/Leaf/007.png
Sigma=4

time stack test :EquivarianceDiagram-test --test-arguments "${InputPath} ${Sigma}"
