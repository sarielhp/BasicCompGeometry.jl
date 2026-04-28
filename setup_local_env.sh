#!/bin/bash

# Configuration: Update this path to the absolute path of your package
PACKAGE_PATH="/home/sariel/prog/26/BasicCompGeometry"

echo "Setting up Julia environment in $(pwd)..."

# Run Julia to activate the local project and develop the package
# 1. --project=.  : Activates the environment in the current directory
# 2. Pkg.develop  : Links the local source code instead of downloading it
julia --project=. -e '
using Pkg;
try
    Pkg.develop(PackageSpec(path="'$PACKAGE_PATH'"));
    Pkg.instantiate();
    println("\nSuccess: BasicCompGeometry is now developed in this environment.");
catch e
    println("\nError: Could not add package at '$PACKAGE_PATH'");
    println(e);
    exit(1);
end
'

if [ $? -eq 0 ]; then
    echo "--------------------------------------------------------"
    echo "To use this environment, run: julia --project=."
    echo "Then in Julia, run: using BasicCompGeometry"
    echo "--------------------------------------------------------"
fi
