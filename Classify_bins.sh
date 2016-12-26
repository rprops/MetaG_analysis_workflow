#/bin/bash

# This script classifies selected bins by phylosift.

# Call this script from where you placed the binning output
# and have stored the bin ID file. The bin ID file is a text file
# with on each new line a bin ID (e.g. 304). The fasta files for the bins are assumed
# to be stored in a subdirectory of the main directory (folder/binFolder).

# created by Ruben Props <ruben.props@ugent.be>

#######################################
##### MAKE PARAMETER CHANGES HERE #####
#######################################

# Make sure you named the binning folder according to the
# parameterization (e.g., k4_L3000 for kmer=4 and length threshold=3000)
k=4
L=3000
folder=k${k}_L${L}
binFolder=fasta-bins
input=bins2classify.txt

#fasta extension
ext=fa

# Number of threads
threads=20

####################################################
##### DO NOT MAKE ANY CHANGES BEYOND THIS LINE #####
#####     unless you know what you're doing    #####
####################################################

cat $input | while read ID
do
   echo "[`date`] Starting with bin ${ID}"
   echo ./${folder}/${binFolder}/${ID}.${ext}
   phylosift all --threads $threads --output ./${folder}/phylosift_bin_${ID} ./${folder}/${binFolder}/${ID}.${ext}
done
