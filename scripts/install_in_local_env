#!/bin/bash

# Determine the directory where this script is located
# Since it lives in scripts/, the package root is one level up.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PACKAGE_PATH="$( cd "$SCRIPT_DIR/.." && pwd )"

echo "Setting up Julia environment in $(pwd)..."
echo "Developing package from: $PACKAGE_PATH"

# Run Julia to activate the local project and develop the package
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
