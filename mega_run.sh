#!/bin/bash
make clean
make 

## for each method
for i in {1..4}
do
    ## for each outer timestep
    for j in {1,0.5,0.1,0.05,0.01}
    do 
        echo run_transfer.sh $i $j
        ## for each numODE to time
        for k in {32,128,512,2048,8192,32768,100000}
        do
            echo time_me.sh $i $j $k
        done
    done

done

