#!/bin/bash
#SBATCH -c 1                  # Number of cores on the same node
#SBATCH --partition=igh       # Specifies the partition of the queue
#SBATCH -o slurmlog/slurmlogmytask_%j.out   # File to which STDOUT will be written, %j will be replaced by the jobID
#SBATCH -e slurmlog/slurmlogmytask_%j.err   # File to which STDERR will be written, %j will be replaced by the jobID.

# Loop through each file in the directory
for f in ../Multi-genome_alignments/hg38-centered/* ; do
    file=$(basename "$f")
    #echo "$file : traitement commencé"

    case $file in

    "chr1.maf" | "chr2.maf" | "chr3.maf" | "chr4.maf" | "chr5.maf" | "chr6.maf" | "chr7.maf" | "chr8.maf" | "chr9.maf" | "chr10.maf" | "chrX.maf")
        echo "filename : $file"
        sbatch --mem=40G ./mafReader.sh "$f" UTR/hg38_3p_UTR_coordinates.tsv AGGTAC TCCATG hg38
        ;;

    "chr11.maf" | "chr12.maf" | "chr13.maf" | "chr14.maf" | "chr15.maf" | "chr16.maf" | "chr17.maf")
        echo "filename : $file"
        sbatch --mem=20G ./mafReader.sh "$f" UTR/hg38_3p_UTR_coordinates.tsv AGGTAC TCCATG hg38
        ;;

    "chr18.maf" | "chr19.maf" | "chr20.maf")
        echo "filename : $file"
        sbatch --mem=10G ./mafReader.sh "$f" UTR/hg38_3p_UTR_coordinates.tsv AGGTAC TCCATG hg38
        ;;
    
    "chr21.maf" | "chr22.maf" | "chrY.maf" | "chrM.maf")
        echo "filename : $file"
        sbatch --mem=5G ./mafReader.sh "$f" UTR/hg38_3p_UTR_coordinates.tsv AGGTAC TCCATG hg38

    esac

done
