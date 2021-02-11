#!/usr/bin/bash
set -e

if [ -z "$1" ]; then make ww3_shel; fi

while [ -n "$1" ]; do
    case "$1" in
        -a) make;;
        -c) make clean; make;;
    esac
    shift
done

cd ../RUN/ 
../SEQ_F90/ww3_grid
#nvprof ../SEQ_F90/ww3_shel 2> prof
NV_ACC_NOTIFY=2 nvprof ../SEQ_F90/ww3_shel 2> prof
#nsys profile --force-overwrite true -o profile --stats=true ../SEQ_F90/ww3_shel 2> nsys_prof
#nvprof --metrics dram_read_throughput,dram_write_throughput ../SEQ_F90/ww3_shel 2> metrics
#NV_ACC_NOTIFY=2 ../SEQ_F90/ww3_shel 2> NV_NOTIFY
#../SEQ_F90/ww3_shel

../SEQ_F90/ww3_outp

vimdiff ww3.68060600.spc out/checking_outp/ww3.68060600.spc

cd ../SEQ_F90
