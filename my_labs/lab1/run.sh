#!/bin/sh
EXE=lab1
RESULTS_DIR=results
NUMBER_X_POINT=100
NUMBER_Y_POINT=100
STEPS=5000
OPT_LEVEL="-O0"
PERF_EVENT="-e branch-misses -e L1-dcache-load-misses -e LLC-load-misses"
PROFILE_RUN=false
function print_input() {
echo "Grid: ${NUMBER_X_POINT}x${NUMBER_Y_POINT} Steps:$STEPS"
}
function plot_picture() {
	gnuplot plot_test_png.p 
}

make OPT_LEVEL=$OPT_LEVEL DEBUG_MODE="-DDEBUG_MODE" rebuild
if [ ! -d $RESULTS_DIR ]; then
	mkdir $RESULTS_DIR
fi	
cd $RESULTS_DIR 
echo "Test run"
print_input
../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
echo "Test picture ploting"
plot_picture


cd ..
make OPT_LEVEL=$OPT_LEVEL PROFILE_OPTIONS="-pg" rebuild
cd results

echo "Performance run"
NUMBER_X_POINT=5000
NUMBER_Y_POINT=5000
STEPS=100
print_input
../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 

if [[ PROFILE_RUN == true ]]; then
	NUMBER_X_POINT=5000
	NUMBER_Y_POINT=5000
	STEPS=100
	echo "Gprof run"
	echo "Grid: ${NUMBER_X_POINT}x${NUMBER_Y_POINT} Steps:$STEPS"
	../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
	gprof ../$EXE > ${EXE}.profile.txt
	gprof ../$EXE | gprof2dot -f prof -n0 -e0 | \
	dot -T png  -o callgraph.png
	gprof -A ../$EXE > ${EXE}.source.txt

	cd ..
	make OPT_LEVEL=$OPT_LEVEL PROFILE_OPTIONS="-g" rebuild
	cd results
	NUMBER_X_POINT=5000
	NUMBER_Y_POINT=5000
	STEPS=100
	echo "Perf run: collecting branch, L1 and LCC misses"
	echo "Grid: ${NUMBER_X_POINT}x${NUMBER_Y_POINT} Steps:$STEPS"
	sudo perf record -g $PERF_EVENT ../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
fi



