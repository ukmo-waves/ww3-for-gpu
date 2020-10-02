#!/usr/bin/python

# change_xygrid.py: A simple python script to change the
# dimensions of the grid for ww3_grid. This edits all the
# required lines in ww3_grid.nml and HOMOGENOUS.depth.

# Note: On Isambard there is difficulties with the LANG 
# environment variable being setup for C by default. Add
# the following lines to .bashrc to avoid this issue.
# >    LANG=en_GB.utf8
# >    export LANG

import numpy as np
import sys, getopt

def main(argv):
    if len(argv) == 0: 
        print('change_xygrid.py <number of x/y points>')
        sys.exit(2)
    elif argv[0] == '-h':
        print('change_xygrid.py <number of x/y points>')
        sys.exit(2)
    else:
        n=int(argv[0]) # Set the number of x/y spaces. 
        SX ="  RECT%%NX           = %i \n" %n
        SY ="  RECT%%NY           = %i \n" %n
        s = "1 "*n + "\n"

        # Make changes to the input file for ww3_grid.nml.
        X = open('HOMOGENEOUS.depth','w+')
        for x in np.linspace(1,n,n+1):
                X.write(s)
        X.close()

        # Create a list of strings identical to the ww3_grid.nml
        # and then exchange two lines within the file.
        with open("ww3_grid.nml", 'r') as fp:
                f = fp.readlines()
                f[170] = SX
                f[171] = SY

        # Write the .nml file with the new lines.
        with open("ww3_grid.nml", 'w') as file:
                file.writelines(f)

if __name__ == "__main__":
	main(sys.argv[1:])
