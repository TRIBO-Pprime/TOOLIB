#!/bin/sh

cd    UTILS ; make all
cd ../QSORT ; make all
cd ../CHOLE ; make all
cd ../ALGEN ; make all
cd ../LEAST ; make all
cd ../MINIM ; make all
cd ../SPLIN ; make all
cd ../TCHEV ; make all
cd ../DIGIS ; make all
cd ../FILES ; make all
cd ../GPLOT ; make all
cd ../FFTW3 ; make all
cd ../INTPL ; make all
cd ../MSOLV ; make all

echo "OK"
read wait

