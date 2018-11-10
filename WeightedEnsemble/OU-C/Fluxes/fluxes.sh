#!/bin/bash

make
rm -f Bins.txt
rm -f WEParams.txt

for runs in {1..10}
do

for Z in 10 15 20 25 30
do

cp ${Z}WEParams.txt WEParams.txt
cp ${Z}Bins.txt Bins.txt
./WE_OU simOutZ${Z}Run${runs}.txt fluxOutZ${Z}Run${runs}.txt errFileZ${Z}Run${runs}.txt 0
rm Bins.txt
rm WEParams.txt
cp simOutZ${Z}Run${runs}.txt Data/simOutZ${Z}Run${runs}.txt
rm simOutZ${Z}Run${runs}.txt
cp fluxOutZ${Z}Run${runs}.txt Data/fluxOutZ${Z}Run${runs}.txt
rm fluxOutZ${Z}Run${runs}.txt
cp errFileZ${Z}Run${runs}.txt Data/errFileZ${Z}Run${runs}.txt
rm errFileZ${Z}Run${runs}.txt

done

done