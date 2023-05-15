#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-user=nicolasmoya2024@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output="blastp.oe"
#SBATCH --job-name="blastp"

wkdir="/projects/b1042/AndersenLab/work/nic/pipeline"


source activate blastx

filename="$(basename $library)"


blastp -query /projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/Plots/N2_library/c_elegans.PRJNA13758.WS279.protein_coding.prot.fa -db $library -out ${filename%.*}.recipro.pb.out -outfmt 6 -num_threads 24

