#!/bin/bash
set -e
for i in */; do
cd $i
grep -c HISEQ dt_int.fasta |  awk -v pwd="${PWD##*/}" '{print pwd, $1}' >> ../sample_sizes.txt
cd -
done
