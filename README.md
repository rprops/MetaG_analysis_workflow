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
Install trimmomatic in your path (http://www.usadellab.org/cms/?page=trimmomatic) to remove adapters. Unzip binary (source may give difficulties) in your local directory:
```
unzip Trimmomatic-Src-0.36.zip
```
and add it to your path (.bash_profile):
```
export PATH=/home/yourusername/Trimmomatic-0.36:$PATH
```

### Quality trimming of reads + removal of adapters
Modules required for quality trimming of reads (adapter contamination/read quality).
```
module load Scythe/0.993b sickle/1.33.6
```
Make sure that you have non-interleaved fastq files of forward and reverse reads. These should have an *R1* tag in their filename.
### Start assembly
Create R1.csv and R2.csv files for megahit
```
ls *1.fastq | head -c -1 | tr '\n' ',' > R1.csv
ls *2.fastq | head -c -1 | tr '\n' ',' > R2.csv
```
Assemble the reads
```
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 36 -o Assembly --presets meta-large > megahit.out&
```
