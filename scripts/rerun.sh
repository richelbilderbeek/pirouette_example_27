#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#SBATCH --partition=gelifes
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=pirex27
#SBATCH --output=example_27.log
#
rm -rf example_27
rm errors.png
time Rscript example_27.R
zip -r pirouette_example_27.zip example_27 example_27.R scripts errors.png

