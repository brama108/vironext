#!/bin/bash

# Function to log messages with timestamps
log() {
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $1"
}

# Function to install new Conda environments
install_conda_envs() {
    log "Installing new Conda environments..."
    for file in ./workflow/envs/*.yml; do
        conda env create -f "$file" >>install.log 2>&1 &
    done
    wait
    log "New Conda environments installed."
}

# Main function
main() {
    while true; do
        echo
        echo
        echo "###################################################################################################"
        echo "###################################################################################################"
        echo "This script will install conda environments required to run the PEI viromics pipeline."
        echo "Please make sure to manually delete the following existing conda environments before proceeding."
        for file in ./workflow/envs/*.yml; do
             echo $(basename "$file" .yml)
        done
        echo "Do you agree to proceed with the installation? [y/n]"

        read varname

        if [ "$varname" = "y" ] || [ "$varname" = "Y" ]; then
            log "Welcome!"
            install_conda_envs
            echo "Installation completed."
            exit 0
        elif [ "$varname" = "n" ] || [ "$varname" = "N" ]; then
            log "Installation aborted."
            exit 0
        else
            echo "Please enter 'y' or 'n'."
        fi
    done
}

# Call the main function
main

