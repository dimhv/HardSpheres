#!/bin/bash

# Define the starting value for 'i'
i=2

# Define the number of times you want to execute
num_executions=58

# Loop to execute the code multiple times with 'i' increasing by 5 each time
for ((j = 0; j <= num_executions; j++)); do
    echo -e "\n Executing iteration $j with argument $i\n"
    ./3d_hs_parallel_main "$i" "$i"
    i=$((i + 5))  # Increment 'i' by 5
    wait
done
