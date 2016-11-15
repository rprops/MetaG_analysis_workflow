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
Look [here] (https://github.com/rprops/DESMAN/wiki/Software-installations) to find out if you need to install anything else for the analysis: 

### Quality trimming of reads

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
### Optional: take random subsample from each sample
This can be required for co-assemblies which are too big. <code>-s</code> sets seed for random sampling.
```
seqtk sample -s 777 *.fastq 5000000 > *.fastq
```
### Start co-assembly
At this point you have multiple softwares to choose from (IDBA_UD/Megahit/...). We choose here for IDBA_UD.

#### IDBA_UD assembly

#### Ray assembly

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
