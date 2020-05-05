Batch=4800
Threads=24
NumPoint=64
ThetaFreq=100
RFreq=100
Delta=1
MaxR=64
STD=50
ScaleFactor=1.0
time stack test :LogFourierSeries-test --test-arguments "${NumPoint} ${ThetaFreq} ${RFreq} ${Delta} ${MaxR} ${STD} ${ScaleFactor} ${Batch} ${Threads} +RTS -N${Threads} -RTS"
