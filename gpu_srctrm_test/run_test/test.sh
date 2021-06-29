#!/usr/bin/bash
set -e

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#                                                                         #
# Shell script created to handle running the WaveWatch III source term    #
# miniapp. This has been adapted for multiple calls test.sh should        #
# be run from a python script (run_test.py) on Ismabard using:            #
#                                                                         #
# qsub -I -q pascalq -l walltime=36000  # or voltaq                       #
#                                                                         #
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 

# Load the nameslist and input files for the executables. 
cd ../RUN/ 

ln -sf inp/ww3_grid.nml .
ln -sf inp/ww3_shel.inp .
ln -sf inp/ww3_outp.inp .

# Run grid preprocessor for WW3
../SEQ_F90/ww3_grid > /dev/null 2>&1

# Run wave model
../SEQ_F90/ww3_shel > /dev/null 2>&1

# Run point output
../SEQ_F90/ww3_outp > /dev/null 2>&1

cd ../run_test

# Alternate options for running ww3_shel, these provide more information for diagnostic and development.

# NVidia profiler only
# nvprof ../SEQ_F90/ww3_shel 2> prof

# NVidia profiler and information for the kernal launches(1) and data transfers(2), or both (3). 
# NV_ACC_NOTIFY=3 nvprof ../SEQ_F90/ww3_shel 2> prof

# Full formed Nvidia profile description. --stats generates summary statistics, --force-overwrite will 
# replace all existing result files with same output file name. 
# nsys profile --force-overwrite true -o profile --stats=true ../SEQ_F90/ww3_shel 2> nsys_prof

# NVidia profiler with device memory read and write throughput. 
# nvprof --metrics dram_read_throughput,dram_write_throughput ../SEQ_F90/ww3_shel 2> metrics

# Cuda debugging option, relevant for issues with Isambard or compiler. 
# cuda-gdb ../SEQ_F90/ww3_shel
