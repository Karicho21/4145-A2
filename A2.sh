#!/bin/bash
#SBATCH --job-name=A2
#SBATCH --partition=Centaurus
#SBATCH --time=01:00:00
#SBATCH --mem=16G



./A2 100 1 10000 100
./A2 2 200 5000000 100
./A2 1000 1 10000 100
