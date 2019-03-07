#!/bin/bash
Input=output/test/Resize/triangle.png
Output=out.png
Size=1024
Threads=1

time stack test :Resize-test --test-arguments "${Input} ${Output} ${Size} ${Threads} +RTS -N${Threads} -s -RTS"
