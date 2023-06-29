#!/bin/bash
#SBATCH -c 1 # Number of cores on the same node
#SBATCH --partition=igh # specifies the partiton of the queue
#SBATCH -o slurmlog/slurmlogmytask_%j.out # File to which STDOUT will be written, %j will be replaced by the jobID
#SBATCH -e slurmlog/slurmlogmytask_%j.err # File to which STDERR will be written, %j will be replaced by the jobID.

################################################################################
# This scripts filters mulputiple alignment files .maf using UTR data          #
#                                                                              #
# If problems encoutered concerning the input or usage of the script           #
# please refer to the help section below or use : ./mafReader -h               #
#                                                                              #
#  version updates and perspectives :-                                         #
#                                    -                                         #
#                                    -                                         #
#                                    -                                         #
#                                    -                                         #
################################################################################

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   RED='\033[1;33m'
   NC='\033[0m'

   echo "####This script takes 5 arguments as follow and in this order"
   echo
   echo
   echo -e "${RED}Usage : ${NC}./mafReader.sh mafFile UTRFile SeedMatchStrand+ SeedMatchStrand- SpecieCentered"
   echo
   echo -e "${RED}-mafFile :${NC} obtained from UCSC genome browser data"
   echo -e "${RED}-UTRFile :${NC} containing all the UTR positions of genes of the specie that the mafFile is centered on"
   echo -e "${RED}-SeedMatchStrand :${NC} Target sequence of the desired miRNA, the seedMatch sequence of both strands is needed"
   echo -e "${RED}-SpecieCentered :${NC} ID of the ${RED}specie${NC}  that the maf file is centered on : for example 'hg38' or 'galGal6'"

}


################################################################################
# Process the input options. Add options as needed.                            #
################################################################################

# Get the options
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
     \?) # incorrect option
         echo "Error: Invalid option"
         exit;;
   esac
done



################################################################################
#                              Main program                                    #
################################################################################


# my arguments #
maf_file=$1
UTR_file=$2
sequence_pos=$3
sequence_neg=$4
species=$5

# my counts #
found_specie=0
found_block=0
found_UTR=0

# split the filename from the Path #
p=$maf_file
file=$( echo ${p##*/} )

# extract the number of the chromosome from the file name : for now I only take in count chrXX.maf and chrUnXXYYY.maf name types #
if [[ "$file" = chr* ]]; then
    #echo "$file name's starts with chr !"
    Current_Chrom="${file#chr}"
    Current_Chrom="${Current_Chrom%.maf}"  
    #echo "$Current_Chrom"
elif [[ "$file" = chrUn* ]]; then
    #echo "$file name's starts with chrUn !"
    Current_Chrom="${file#chrUn}"
    Current_Chrom="${Current_Chrom%.maf}"  
    #echo "$Current_Chrom"
#else
    #echo "$file name's does not match what I expect !! "
fi


# Iterate over each line in the maf file #
echo "Le fichier traité ici est $file"
while IFS= read -r line
do
    # Check if the line starts with 'a' indicating a NEW alignment block #
    if [[ $line == "a"* ]]; then
        found_block=0
        found_specie=0
        found_UTR=0
    fi
    # Check if the line starts with 's' indicating a sequence line
    if [[ $line == "s"* ]]; then
        # Extract the sequence from the line using an array list #
        #arrayTest=("tata" "toto" "titi")
        #echo ${arrayTest[1]}
        line_parts=($line)

        #this variable contains for each line the sequence that is aligned
        current_sequence=${line_parts[6]//-/}

        #this variable contains for each line starting position of the sequence
        position=${line_parts[2]}

        #this variable contains for each line the name of the specie
        source_genome=${line_parts[1]}

        shopt -s nocasematch;
        # Check if the given specie is found in the current alignment block #
        if [[ $source_genome == *"$species"* ]]; then
        found_specie=1

            # Check if in this line where the specie was found we have at least ONE seedmatch (+ or - or both) #
            
            if [[ $current_sequence == *"$sequence_pos"* ]] || [[ $current_sequence == *"$sequence_neg"* ]]; then
            found_block=1
            #echo "I found a block with the good specie and at least on seed match (+ or -)"

            # If the two precedent conditions are met, iterate through the UTR file to see if the sequence is present in an UTR region found on that same chromosome #
            while IFS= read -r line2
                do
                    line_parts2=($line2)

                    # This variable contains the chromosome number/ID in the UTR file, so it can be compared with the chromosome number(Current_Chrom) of interest
                    chrom=${line_parts2[3]}

                    #Start position of the gene in the UTR
                    start_pos=${line_parts2[5]}

                    #End position of the gene in the UTR
                    end_pos=${line_parts2[6]}

                    #Check if the line in the UTR file is about the chromosome that we are searching for
                    if [[ $chrom == $Current_Chrom ]]; then
                                        #echo "I am in the right section of the UTR from chrom : $chrom AND $Current_Chrom, I search for matches"

                        #If it is the needed chromosome number, we check if the start position of the sequence in the MAF file is between the gene positions.
                        if [[ $start_pos -le $position && $end_pos -ge $position ]]; then
                            found_UTR=1
                            echo  "###########################################################################################"    
                        fi 
                    fi    
                done < "$UTR_file"

                
            fi
        fi
    fi

    # Print the line if the sequence was found in the alignment block and corresponds to the specie  and corresponds to a valid UTR #
    if [[ $found_block -eq 1 ]] && [[ $found_specie -eq 1 ]] && [[ $found_UTR -eq 1 ]] ; then
        if [[ $line == "s"* ]]; then
            echo "$line" 
        fi
    fi
done < "$maf_file"

