#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
sys.path.insert(0,'/home/mo-lsampson/ww3-for-gpu/gpu_srctrm_test/RUN/inp')
import change_xygrid
import time

def main():
    # The main purpose of this routine is to produce a graphic for the time 
    # comparisons between different model runs of either, ACC, OMP or SEQ. 
    # This will run for each grid size (N^2) and each modelruns variable. 
    gridsizes = [5, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 240, 280, 320]
    rundir = '/home/mo-lsampson/ww3-for-gpu/gpu_srctrm_test/RUN'
    pwd=os.getcwd()

    # We are only doing MM for the ACC run. Hardcoded in makefile. It is worth being
    # aware that the SEQ runs are much longer than the other 2. 
    modelruns = []#'Sequential','ACC run','OMP run'] # 'Sequential','ACC run','OMP run'

    # The script will attempt to save the current time outputs as much as possible
    # so there is no need to re-run completely if a failure occurs. Instead remove the
    # successful runs from the modelruns list. Requires consistent gridsizes. 

    colors=['k','g','r']

    duration = np.zeros((len(gridsizes), 3))
    timer, seq, acc, omp = np.zeros((len(gridsizes))),np.zeros((len(gridsizes))),\
            np.zeros((len(gridsizes))),np.zeros((len(gridsizes)))
    for j, mr in enumerate(modelruns):
        # For the OMP run we use the OMP_F90 directory. This has not been changed
        # since the github was created, other than to allow for testing. 
        timer[:]=0.0
        if mr == 'OMP run':
            srcdir='/home/mo-lsampson/ww3-for-gpu/gpu_srctrm_test/OMP_F90'
        else:
            srcdir='/home/mo-lsampson/ww3-for-gpu/gpu_srctrm_test/SEQ_F90'

        # We clean the make file to ensure since both SEQ and ACC will make from
        # the same list of .F90 files in same directory.
        os.chdir(srcdir)
        clean = subprocess.Popen(["make clean"], shell=True)
        clean.wait()
        
        # make -s holts any output from the make file. This is not useful for 
        # development but does help here. MODE is defined in the make file to
        # select the right FFLAGS that are required. 
        mkstr="make -s MODE="+mr[0:3].lower()
        process=subprocess.Popen([mkstr], shell=True)
        process.wait() # Without wait the python script runs on to next tasks.


        # To keep all runs consistent in terms of runtime comparison we premake
        # the model definitions from a clean slate. This means the first call
        # to test_shel.sh will skip the ww3_grid (or will be very fast). 
        os.chdir(rundir)
        rmst = subprocess.Popen(["rm ST4TABUHF2.bin"], shell=True)
        rmst.wait()
        os.system(srcdir+"/ww3_grid")
        os.chdir(pwd)

        for i, gs in enumerate(gridsizes):
             print(gs)
             os.chdir(rundir+'/inp')
             change_xygrid.main([str(gs)])
             os.chdir(pwd)
             s = time.time()
             subprocess.call(['sh','./test.sh'])
             e = time.time()
             timer[i]=(e-s)
        if mr == 'OMP run':
            omp=timer.copy()
        elif mr == 'Sequential':
            seq=timer.copy()
        elif mr == 'ACC run':
            acc=timer.copy()
        else:
            print('Error in pre-assigned model run names')
            exit(0)

    try:
        temp = np.load('../run_test/timecost.npz')
        duration[:,0] = temp['seq']
        duration[:,1] = temp['acc']
        duration[:,2] = temp['omp']
        print('Loaded timecost.npz')
        temp.close()
    except:
        print('Error in loading timecost.npz. Saving currently available arrays.')
        np.savez('../run_test/timecost', seq=seq, acc=acc, omp=omp)

    for k, mode in enumerate((seq, acc, omp)):
        if all(mode == np.zeros((len(gridsizes)))):
            print(mode)
            mode = duration[:,k]
            print(mode)
        else:
            print(mode)
            duration[:,k]=mode

    np.savez('../run_test/timecost.npz', seq=duration[:,0], acc=duration[:,1], omp=duration[:,2])
    #plt.title('Runtime for the mini-app using different compilation options')
    plt.xlabel('Number of gridpoints (N^2)')
    plt.ylabel('Time (Seconds)')
    for j, mr in enumerate(('Sequential','ACC run','OMP run')):
        plt.plot(gridsizes,duration[:,j],color=colors[j],marker='^',linestyle='--',label=mr)
    plt.legend()
    plt.savefig('../run_test/timecomp.png')
    #plt.show()


if __name__ == '__main__':
    main()
