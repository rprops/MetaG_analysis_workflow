# DESMAN

```R
module load gsl python-anaconda2
cd
git clone https://github.com/chrisquince/DESMAN.git
cd DESMAN
CFLAGS="-I$GSL_INCLUDE -L$GSL_LIB" python setup.py install --user
```
It is now installed and can be called like:
```R
~/.local/bin/desman -h
```
Or add it to your path:
```R
export PATH=~/.local/bin:$PATH
```
