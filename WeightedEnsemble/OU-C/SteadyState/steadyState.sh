#!/bin/bash

make
for runs in {1..10}
do
./WE_OU simOut${runs}.txt fluxOut${runs}.txt errFile${runs}.txt 0
done