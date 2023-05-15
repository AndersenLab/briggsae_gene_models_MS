#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output=bedcov.oe
#SBATCH --job-name="bedcov"

source activate bedtools

bedtools coverage -sorted -a exon_pos.bed -b ../../../gene_predictions/QX1410/alignments/Aligned.sortedByCoord.out.bam
