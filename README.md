# DESMAN

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

Create R1.csv and R2.csv files for megahit
```
ls *1.fastq | head -c -1 | tr '\n' ',' > R1.csv
ls *2.fastq | head -c -1 | tr '\n' ',' > R2.csv
```
