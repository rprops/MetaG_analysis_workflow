# DESMAN

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
Look [here] (https://github.com/rprops/DESMAN/wiki/Software-installations) to find out if you need to install anything else for the analysis. Make sure all these modules are loaded.

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
At this point you have multiple softwares to choose from (IDBA-UD/Megahit/...). We choose here for IDBA-UD.

#### IDBA-UD assembly
Merge all interleaved files to one file with the following shell script:
```
bash assembly_prep.sh
```
**Optional:** normalize data based on coverage (BBNorm)
This can be required for co-assemblies which are too big (see: [here] (http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide/)). First estimate memory requirements based on number of unique kmers using <code>loglog.sh</code>. Run this in a pbs script (due to memory requirements).
```
loglog.sh in=merged_dt_int.fasta
bbnorm.sh in=merged_dt_int.fasta out=merged_dt_int_normalized.fasta target=100 min=5
```
Due to the normalization it can happen that some paired reads may have lost one of their mates. Therefore we must split the normalized fasta file and recombine them to a single fasta file (without the unpaired sequences):
```
# Script courtesy of https://github.com/MadsAlbertsen/miscperlscripts
# Split interleaved fasta
/nfs/vdenef-lab/Shared/Ruben/scripts_metaG/SeqTools/splitpe.fasta.pl merged_dt_int_normalized.100.5.fasta
# Interleave fasta
/nfs/vdenef-lab/Shared/Ruben/scripts_metaG/SeqTools/interleave.pl -fwd p1.fa -rev p2.fa -o merged_dt_int_normalized.100.5_pe
```
Now run idba_ud
```

```



#### Ray assembly
Make sure <code>gcc</code> and <code>openmpi</code> modules are loaded.

#### Megahit assembly
Put all your fastq files from all samples in one folder. Create R1.csv and R2.csv files for Megahit:
```
ls *.R1.fastq | head -c -1 | tr '\n' ',' > R1.csv
ls *.R2.fastq | head -c -1 | tr '\n' ',' > R2.csv
```
Assemble the reads (do not run this on login nodes).
```
megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 36 -o Assembly --presets meta-large > megahit.out
```
