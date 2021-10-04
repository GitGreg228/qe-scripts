#!/bin/sh

for dname in $(ls -d */)
do
cd $dname
sbatch script1.sh
cd ..
done
