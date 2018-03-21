#!/bin/bash
#make clean
#make 
./POSTERR sample_input.txt output/$1_$2_output.txt $1 $2
scp output/$1_$2_output.txt abg6257@quest.northwestern.edu:cudaCODE/output
