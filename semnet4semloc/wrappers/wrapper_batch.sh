#!/bin/bash

echo "PBS_ARRAYID " $PBS_ARRAYID
echo $myvar

# myvar: python script to execute
# PBS_ARRAYID: subject index to be processed

/imaging/local/software/anaconda/2.4.1/2/bin/python2.7 $myvar $PBS_ARRAYID
