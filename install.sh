#!/bin/bash
cd "${0%/*}" || exit                                # Run from this directory
#------------------------------------------------------------------------------

# Define the default symlink location
SYMLINK_PATH="$WM_PROJECT_USER_DIR/src/simlinks"
LIB_USER_SRC="$PWD/src"

# Check if the src directory exists
if [[ ! -d "$LIB_USER_SRC" ]]; then
    echo "Error: Directory '$LIB_USER_SRC' does not exist. Please check the path and try again."
    exit 1
fi

# Create the symlink directory if it doesn't exist
mkdir -p "$SYMLINK_PATH"

# Search for all lnInclude directories under LIB_USER_SRC
find "$LIB_USER_SRC" -type d -name "lnInclude" | while read -r lnIncludeDir; do
    # Extract the library name
    libName=$(basename "$(dirname "$lnIncludeDir")")

    # Create a symlink to this lnInclude directory
    ln -sf "$lnIncludeDir" "$SYMLINK_PATH/$libName"

    # Print created OpenFOAM library symlinks
    echo "Linked $lnIncludeDir -> $SYMLINK_PATH/$libName"
done

# Compile user libraries and solvers
chmod u+x src/Allwmake
chmod u+x applications/solvers/Allwmake
bash src/Allwmake
bash applications/solvers/Allwmake
