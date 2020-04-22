#!/bin/sh
EXE=lab1
RESULTS_DIR=results
NUMBER_X_POINT=600
NUMBER_Y_POINT=600
STEPS=1500
make rebuild
if [ ! -d $RESULTS_DIR ]; then
	mkdir $RESULTS_DIR
fi	
cd $RESULTS_DIR 
#Performance test run
time ../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
gnuplot plot_test_png.p 

cd ..
make PROFILE_OPTIONS="-g -pg" rebuild
cd results
../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 

gprof ../$EXE > ${EXE}.profile.txt
gprof ../$EXE | gprof2dot -f prof -n0 -e0 | \
dot -T png  -o callgraph.png
gprof -A ../$EXE > ${EXE}.source.txt


