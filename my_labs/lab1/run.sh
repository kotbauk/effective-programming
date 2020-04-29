#!/bin/bash
EXE=lab1
RESULTS_DIR=results
NUMBER_X_POINT=100
NUMBER_Y_POINT=100
STEPS=5000
OPT_LEVEL="-O3 -march=native"
PERF_EVENT="branch-misses,L1-dcache-load-misses,LLC-load-misses,cycles,instructions,cache-misses"
PROFILE_RUN=false
function print_input() {
echo "Grid: ${NUMBER_X_POINT}x${NUMBER_Y_POINT} Steps:$STEPS"
}
function plot_picture() {
	gnuplot plot_test_png.p 
}

make OPT_LEVEL="$OPT_LEVEL" DEBUG_MODE="-DDEBUG_MODE" rebuild
if [ ! -d $RESULTS_DIR ]; then
	mkdir $RESULTS_DIR
fi	
cd $RESULTS_DIR 
echo "Test run"
print_input
../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
echo "Test picture ploting"
plot_picture

cd -
make OPT_LEVEL="$OPT_LEVEL" rebuild
cd $RESULTS_DIR

echo "Performance run"
NUMBER_X_POINT=5000
NUMBER_Y_POINT=5000
STEPS=100
print_input
time ../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
PROFILE_RUN="true"
if [[ $PROFILE_RUN == true ]]; then
	NUMBER_X_POINT=5000
	NUMBER_Y_POINT=5000
	STEPS=100
	cd -
	make OPT_LEVEL="$OPT_LEVEL" PROFILE_OPTIONS="-g -pg" rebuild
	cd $RESULTS_DIR
	echo "Gprof run"
	echo "Grid: ${NUMBER_X_POINT}x${NUMBER_Y_POINT} Steps:$STEPS"
	../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
	gprof ../$EXE > ${EXE}.profile.txt
	gprof ../$EXE | gprof2dot -f prof -n0 -e0 | \
	dot -T png  -o callgraph.png
	gprof -A ../$EXE > ${EXE}.source.txt

	cd -
	make OPT_LEVEL="$OPT_LEVEL" PROFILE_OPTIONS="-ggdb" rebuild
	cd $RESULTS_DIR
	NUMBER_X_POINT=5000
	NUMBER_Y_POINT=5000
	STEPS=100
	echo "Perf run: collecting $PERF_EVENT"
	echo "Grid: ${NUMBER_X_POINT}x${NUMBER_Y_POINT} Steps:$STEPS"
	perf record -e $PERF_EVENT ../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS 
	perf stat -d -d -d -B -e $PERF_EVENT ../$EXE $NUMBER_X_POINT $NUMBER_Y_POINT $STEPS
fi


