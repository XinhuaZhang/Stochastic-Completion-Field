#!/bin/bash
ImageName=triangle
Input=output/test/Resize/${ImageName}.png #Ehrenstein_only.jpeg #triangle.png
Size=128
Output=${ImageName}_${Size}.png
Threads=1

time stack test :Resize-test --test-arguments "${Input} ${Output} ${Size} ${Threads} +RTS -N${Threads} -s -RTS"
