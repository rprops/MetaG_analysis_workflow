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

### General software installations
Modules required for quality trimming of reads (adapter contamination/read quality). 

Install trimmomatic in your path (http://www.usadellab.org/cms/?page=trimmomatic) to remove adapters. Unzip binary (source may give difficulties) in your local directory:
```
unzip Trimmomatic-Src-0.36.zip
```
and add it to your path (.bash_profile):
```
export PATH=/home/yourusername/Trimmomatic-0.36:$PATH
```
Check if permissions allow execution.
```
cd ./Trimmomatic-0.36
ls -l
```
```
module load Scythe/0.993b sickle/1.33.6
```
Scythe is not installed on flux so install it in your home directory and add it to the path similar to Trimmomatic. 
```
cd
git clone https://github.com/vsbuffalo/scythe.git
cd scythe
make all
```
Also install Fastqc binary from here: http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc. And add to path (see trimmomatic) and set permissions.
```
cd ./FastQC
chmod 755 fastqc
```
Make sure to reconnect to the server in order to update your path.

### Quality trimming of reads

Make sure that you have non-interleaved fastq.gz files of forward and reverse reads. These should have an *R1* tag in their filename and saved in a directory called *sample*, e.g.: sample/sample.R1.fastq.gz.

**IMPORTANT** The adapter trimming in qc.sh is standard set for Truseq paired-end libraries. You will have to adjust this if this is different for you by changing the path in the shell script to the correct fasta file of the adapters located at /home/your_username/Trimmomatic-0.36/adapters.
```
bash /nfs/vdenef-lab/Shared/Ruben/scripts_metaG/wrappers/Assembly/qc.sh sample_directory
```
Alternatively modify the <code>run_quality.sh</code> or <code>run_quality.pbs</code> scripts to run sequentially.

### Start assembly
At this point you have multiple softwares to choose from (IDBA_UD/SPAdes/Megahit/...). We choose here for IDBA_UD.

#### IDBA_UD assembly

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
