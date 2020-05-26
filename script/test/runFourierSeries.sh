GPU=1
Size=63   #127
Delta=0.125
Period=63 #127
Freq=31   #127
AngularFreq=5
RadialFreq=5
STD=16     #63
#+RTS -N1 -RTS
time stack test :FourierSeries-test --test-arguments "${GPU} ${Size} ${Delta} ${Period} ${Freq} ${AngularFreq} ${RadialFreq} ${STD} +RTS -N28 -RTS "
