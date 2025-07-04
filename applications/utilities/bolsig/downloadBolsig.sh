#!/bin/bash
# ------------------------------------------------------------------------------
# --- README
# Run this bash file to download the 07/2024 version of the BOLSIG+ executable
# for Linux-based operating systems. You only need to do this once. This bash
# file must always be located at the project's root directory.
# 
# --- TROUBLESHOOTING
# The unzip command is usually preinstalled in most Linux distributions.
# If not, then for Ubuntu or Debian-based distributions execute:
#   $ sudo apt update
#   $ sudo apt -y install unzip
# ------------------------------------------------------------------------------
# Run from this directory
cd "${0%/*}" || exit

# Define the download URL of the BOLSIG+ executable
bolsigLinkAdress=https://www.bolsig.laplace.univ-tlse.fr/wp-content/uploads/2024/07/bolsigplus072024-linux.zip

# Download and extract the files
mkdir -p ./runBolsig/bolsigExe
wget -O ./runBolsig/bolsigExe/bolsigFiles.zip $bolsigLinkAdress
unzip -j -o ./runBolsig/bolsigExe/bolsigFiles.zip bolsigminus -d ./runBolsig/bolsigExe
chmod +x ./runBolsig/bolsigExe/bolsigminus
rm ./runBolsig/bolsigExe/bolsigFiles.zip

# Print final message
echo
echo BOLSIG+ executable was downloaded successfully!
