#!/bin/bash
Input=output/test/Resize/Ehrenstein_only.jpeg #Ehrenstein_only.jpeg #triangle.png
Output=out.png
Size=128
Threads=1

time stack test :Resize-test --test-arguments "${Input} ${Output} ${Size} ${Threads} +RTS -N${Threads} -s -RTS"
