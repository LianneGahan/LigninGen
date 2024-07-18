#!/bin/bash
#$1 : number of cores to be used
#$2 : Number of runs to use in a sample
#$3 : CHOICE: 0 to carry out a new gridsearch, 1 to use pre-existing gridsearch data, 
#			2 for the user to choose the starting point
#$4 : Number of samples per degree of freedom in the gridsearch** ignored if $3 == 0
    
if [ "$3" -eq 0 ]; then

	# Perform a grid search in the first instance to choose the starting conditions
	mkdir Gridsearch/
	(echo "Making folders for grid search")
	# Sample the parameter space to choose starting conditions
	N_varied_params=$(grep -o '1' kin_params_to_randomize.txt | wc -l)
	echo "------------------------------------"
	echo "Conducting a grid search with: "
	echo "Number of varied parameters: $N_varied_params"
	Nsamples=$(($4**($N_varied_params)))
	echo "Total number of samples: $Nsamples"
	echo "------------------------------------"
	if [ "$4" -eq 0 ]; then
		Nsamples=1
	fi
	for j in $(seq 1 $Nsamples)
	do
		mkdir Gridsearch/Run_${j}
		cp -a latest/Output/ Gridsearch/Run_${j}/
	done
	
	(echo "Changing initial conditions")
	(echo "Running")

	python3 latest/Code/grid_search.py $4 $1 #>/dev/null & 
	echo "Grid Search Complete"
	PROCESSID=$!
	wait $PROCESSID
fi

if [ "$3" -eq 1 ]; then
	# Run python script to read in parameters
	python3 readin_gridsearch.py  #>/dev/null &
	PROCESSID=$!
	wait $PROCESSID
	# Need to make a function which converts a grid search into this format
fi
