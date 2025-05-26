#!/bin/bash

# Define arrays
OBSID=('ASPeCt' 'SB_AWI')
# AEM_AWI
SATELLITE=('CS2' 'ENV')
HS=('NH' 'SH')

# Max parallel jobs
MAX_JOBS=10

# Function to limit the number of parallel jobs
wait_for_jobs() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 1
    done
}

for id in "${OBSID[@]}"; do
    script_dir="/dmidata/users/ilo/projects/RRDPp/code/$id"

    if [ -d "$script_dir" ]; then
        cd "$script_dir" || continue

        # Find all python scripts that match *processing_final.py
        scripts=$(find . -type f -name '*processing_final.py')

        if [ -n "$scripts" ]; then
            for script in $scripts; do
                wait_for_jobs  # Check before launching new job
                echo "Running $script for $id..."
                python -W ignore "$script" &
            done
        else
            echo "No *processing_final.py scripts found in $script_dir"
        fi
    else
        echo "Directory $script_dir does not exist"
    fi
done

# Wait for all background jobs to finish
wait

