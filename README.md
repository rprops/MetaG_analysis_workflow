# Metagenomic analysis workflow

### Install DESMAN
```
module load gsl python-anaconda2
cd
git clone https://github.com/chrisquince/DESMAN.git
cd DESMAN
CFLAGS="-I$GSL_INCLUDE -L$GSL_LIB" python setup.py install --user
```
It is now installed and can be called like:
```
~/.local/bin/desman -h
```
Or add it to your path:
```
export PATH=~/.local/bin:$PATH
```
### Other software
Look [here](https://github.com/rprops/DESMAN/wiki/Software-installations) to find out if you need to install anything else for the analysis. Make sure all these modules are loaded.

### Step 1: Quality trimming of reads

Make sure that you have non-interleaved fastq.gz files of forward and reverse reads. These should have an *R1* tag in their filename and saved in a directory called *sample*, e.g.: sample/sample.R1.fastq.gz.

**IMPORTANT** The adapter trimming in qc.sh is standard set for Truseq paired-end libraries. You will have to adjust this if this is different for you by changing the path in the shell script to the correct fasta file of the adapters located at /home/your_username/Trimmomatic-0.36/adapters.

**IMPORTANT** In case you are unsure which adapters are present in the sequences, you can download bbtools
(https://sourceforge.net/projects/bbmap/) unzip the tar.gz and add the directory to your export path. Then run the following code on a subsample of a sample (e.g., 1m reads). The resulting consensus sequences of the adapters will be stored in <code>adapters.fa</code>. 
```
bbmerge.sh in1=$(echo *R1.fastq) in2=$(echo *R2.fastq) outa=adapters.fa reads=1m
```
Run quality trimming:
```
bash /nfs/vdenef-lab/Shared/Ruben/scripts_metaG/wrappers/Assembly/qc.sh sample_directory
```
Alternatively modify the <code>run_quality.sh</code> or <code>run_quality.pbs</code> scripts to run sequentially.

Copy Fastqc files to new folder (adjust the paths)
```
rsync -a --include '*/' --include '*fastqc.html' --exclude '*' /scratch/vdenef_fluxm/rprops/DESMAN/metaG/Nextera /scratch/vdenef_fluxm/rprops/DESMAN/metaG/FASTQC --progress
```
#### Random subsampling of reads. Works if you are interested in abundant taxa.

You can check the number of reads in the interleaved fasta file with the <code>sample_size.sh</code>. This will store the sample sizes in the <code>sample_sizes.txt</code> file. Run this in the directory where your samples are located. Then use seqtk to randomly subsample your fasta files to the desired sample size. <code>-s</code> sets seed for random sampling.
```
bash sample_size.sh
seqtk sample -s 777 *.fastq 5000000 > *.fastq
```
### Step 2: Start co-assembly
At this point you have multiple assemblers to choose from (IDBA-UD/Megahit/...). We choose here for IDBA-UD.

#### IDBA-UD assembly
Merge all interleaved files to one file with the following shell script:
```
bash assembly_prep.sh
```
**Optional:** normalize data based on coverage (BBNorm)
This can be required for co-assemblies which are too big (see: [here] (http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide/)). First estimate memory requirements based on number of unique kmers using <code>loglog.sh</code>. Run this in a pbs script (due to memory requirements).
```
loglog.sh in=merged_dt_int.fasta
bbnorm.sh in=merged_dt_int.fasta out=merged_dt_int_normalized.100.5.fasta target=100 min=5
```
Due to the normalization some paired reads will have lost their mate. Therefore we must split the normalized fasta file and recombine them to a single fasta file (without the unpaired sequences). Run the below code in a pbs script. Adjust the pattern <code>:1</code> or <code>:2</code> to reflect the R1/R2 read names in your sequence data. Also adjust 
```
# Generate list of R1 and R2 reads
grep " 1:" merged_dt_int_normalized.100.5.fasta | sed "s/>//g" > list.R1
grep " 2:" merged_dt_int_normalized.100.5.fasta | sed "s/>//g" > list.R2
# Extract R1 and R2 reads from interleaved file using the filterbyname.sh script from BBmap
filterbyname.sh in=merged_dt_int_normalized.100.5.fasta names=list.R1 out=merged_dt_int_normalized.100.5.R1.fasta -include t
filterbyname.sh in=merged_dt_int_normalized.100.5.fasta names=list.R2 out=merged_dt_int_normalized.100.5.R2.fasta -include t
# Interleave R1/R2 reads using inhouse perl script (author: Sunit Jain)
/nfs/vdenef-lab/Shared/Ruben/scripts_metaG/SeqTools/interleave.pl -fwd merged_dt_int_normalized.100.5.R1.fasta -rev merged_dt_int_normalized.100.5.R2.fasta -o remerged_dt_int_normalized.100.5
```

Now run idba_ud (you may have to adjust the parameters.
```
idba_ud -o idba_k52_100_s8 -r remerged_dt_int_normalized.100.5.fasta --num_threads ${PBS_NP} --mink 52 --maxk 100 --step 8
```
In case you run this on flux you'll have to add the following line of code above your actual code to allow openmpi to run multithreaded.
```
# On Flux, to use OpenMP, you need to explicitly tell OpenMP how many threads to use.
export OMP_NUM_THREADS=${PBS_NP}
```

### IMPORTANT
Some assemblers (such as the newest idba version) will give extra information in the headers which can be problematic for other software (e.g anvio and phylosift). To avoid this we can use an anvio script to make the headers compatible with all software.
```
module load anvio
anvi-script-reformat-fasta contig.fa -o contigs-fixed.fa -l 0 --simplify-names
mv contigs-fixed.fa contig.fa
```

### Step 3: Generating assembly stats
We use quast for this purpose (adjust path for your specific assembly):
```
quast.py -f --meta -t 20 -l "Contigs"  megahit_assembly_sensitive/final.contigs.fa
```
For IDBA-UD this can be run both on the scaffolds and the contigs. In general, the contigs are more reliable.

### Step 4A: CONCOCT binning
CONCOCT binning is not optimal for large co-assemblies since it combines both coverage and kmer/GC% in one step. This will result in very slow clustering for complex data sets. So go for Binsanity, if this is the case. We're gonna focus on the IDBA-UD assembly now but for megahit you can find the analogy [here](https://github.com/rprops/DESMAN/wiki/Megahit-analysis).

#### IDBA-UD
IDBA-UD will also make scaffolds based on the contigs. So the two primary output files are:
```
scaffold.fa
contig.fa
```
First cut your large contigs before running CONCOCT (make sure CONCOCT is in your path).
```
module load python-anaconda2/201607
module load gsl

mkdir contigs
cut_up_fasta.py -c 10000 -o 0 -m idba_k52_100_s8/contig.fa > contigs/final_contigs_c10K.fa
```
Be aware that this introduces dots in the contig names (e.g., contig-100_1647.1). This can pose an issue for taxonomic classifications through Phylosift. Replace dots by underscores as follows:
```
sed "s/\./_/g" final_contigs_c10K.fa > final_contigs_c10K_nodots.fa
```
**Optional/Necessary** You can at this point already remove short contigs (e.g., < 1000) this will save time in the mapping but also for generating the coverage file.
```
reformat.sh in=final_contigs_c10K.fa out=final_contigs_c10K_1000.fa minlength=1000
```
Make an index for your contigs fasta.
```
cd contigs
bwa index final_contigs_c10K.fa
cd -
```
Then map reads (**warning:** fastq files must not contain unpaired reads).
```
for file in *R1.fastq
do

   stub=${file%_R1.fastq}

   echo $stub

   file2=${stub}_R2.fastq

   bwa mem -t 20 contigs/final_contigs_c10K.fa $file $file2 > Map/${stub}.sam
done

```
After the mapping create .bam files using samtools. Make sure you install bedtools2 for the next step.
**IMPORTANT:**Make sure you have samtools >v1.3!
Use <code>-@</code> to multithread and <code>-m</code> to adjust memory usage.
```
#/bin/bash

set -e

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub
    samtools view -h -b -S $file > ${stub}.bam
    samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam
    samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam -@ 8 
    samtools index ${stub}.mapped.sorted.bam
    bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g contigs/final_contigs_c10K.len > ${stub}_cov.txt
done

for i in Map/*_cov.txt 
do 
   echo $i
   stub=${i%_cov.txt}
   stub=${stub#Map\/}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > Map/${stub}_cov.csv
done

~/DESMAN/scripts/Collate.pl Map | tr "," "\t" > Coverage.tsv
```
Once you've formatted the coverage files we can start binning using CONCOCT (use 40 core in pbs script to maximize on C-CONCOCT). We perform the binning on contigs>3000bp and put the kmer-signature on 4.
```
module load gsl
mkdir Concoct
cd Concoct
mv ../Coverage.tsv .
concoct --coverage_file Coverage.tsv --composition_file ../contigs/final_contigs_c10K.fa -s 777 --no_original_data -l 3000 -k 4
cd ..
```
After binning we can quickly evaluate the clusters:
```
mkdir evaluation-output
Rscript ~/CONCOCT/scripts/ClusterPlot.R -c clustering_gt3000.csv -p PCA_transformed_data_gt3000.csv -m pca_means_gt3000.csv -r pca_variances_gt3000_dim -l -o evaluation-output/ClusterPlot.pdf
```
CONCOCT only outputs a list of bins with the associated contig ids. We further use the functionality of the Phylosift software to extract the corresponding fasta file for each bin (run in pbs script). We needs these fasta files later on for CheckM.
```
mkdir fasta-bins
extract_fasta_bins.py ../contigs/final_contigs_c10K.fa ./k4_L3000_diginorm/clustering_gt2000.csv --output_path ./k4_L3000_diginorm/evaluation-output
```
Then run CheckM, this requires pplacer, hmmer and prodigal to be loaded. This will process wil generate an output file <code>bin_stats_mega_k4_L3000.tsv</code> and a plot describing these results, which will be stored in the evaluation-output folder.
```
checkm lineage_wf --pplacer_threads 10 -t 10 -x fa -f bin_stats_mega_k4_L3000.tsv ./k4_L3000/fasta-bins ./k4_L3000/tree_folder
checkm bin_qa_plot ./k4_L3000/tree_folder ./k4_L3000/fasta-bins ./k4_L3000/evaluation-output -x fa --dpi 250
```
Next we select some bins and classify them with the following bash script (classify_bins.sh):
```
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
```
