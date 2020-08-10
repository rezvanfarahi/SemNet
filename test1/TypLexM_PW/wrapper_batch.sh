#!/bin/bash

echo "PBS_ARRAYID " $PBS_ARRAYID
echo $myvar

# myvar: python script to execute
# PBS_ARRAYID: subject index to be processed

/imaging/local/software/EPD/latest/x86_64/bin/python $myvar $PBS_ARRAYID
