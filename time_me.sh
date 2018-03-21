#!/bin/bash
#make clean
#make 
./POSTERR sample_input.txt timings/$3_timings/$1_$2_output.txt $1 $2 $3
scp timings/$3_timings/$1_$2_output.txt abg6257@quest.northwestern.edu:cudaCODE/timings/$3_timings
