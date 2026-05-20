#!/bin/bash

RESULTS_DIR=""

if then
monet --input "Interactome_weighted.txt" --method M1 --largest 80 --output "${RESULTS_DIR}/Unchanged_modules/Interactome_weighted_modules.txt"

monet --input "perturbed_interactome.txt" --method M1 --largest 80 --output "${RESULTS_DIR}/Perturbed_modules/perturbed_interactome_modules.txt"
fi 


for i in 1 2 3 4 5; do
    echo "Running M1 on random_${i}_interactome..."
    monet --input "random_${i}_interactome.txt" --method M1 --largest 80 --output "${RESULTS_DIR}/Random_modules/random_${i}_modules.txt"
done

echo "All networks processed."