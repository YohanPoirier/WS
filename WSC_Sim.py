from __future__ import division
from os import getcwd,path,system,listdir
from fnmatch import filter
from subprocess import Popen
from shutil import copyfile

# The goal of this script is to run several WSC simulations to avoid the memory leak problem.

def Read_State(cwd):

    # cwd    : Current working directory where this subroutine is called.
    # t_Even : Time of writing of State_Even.dat.
    # t_Odd  : Time of writing of State_Odd.dat.

    # This subroutine reads the time in both State_Even.dat and State_odd.dat.

    # Searching State_Even.dat.
    State_Even_matches = []
    files = [f for f in listdir('.') if path.isfile(f)] # Loop over the files of the cwd.
    for file in filter(files, 'State_Even.dat'): # Only *.geom file.
        State_Even_matches.append(path.join(file))
    if(len(State_Even_matches) >= 2):
        print "WSC_Sim: there is more than one State_Even.dat input file. there must be only one!"
    elif(len(State_Even_matches) == 0):
        t_Even = -1
    else: # len(State_Even_matches) = 1
        State_Even = open(State_Even_matches[0],'r')
        State_Even.readline().rstrip() # Title.
        State_Even.readline().rstrip() # Starting time.
        ligne = State_Even.readline().split() # t = time.
        t_Even = float(ligne[0])

    # Searching State_Odd.dat.
    State_Odd_matches = []
    files = [f for f in listdir('.') if path.isfile(f)] # Loop over the files of the cwd.
    for file in filter(files, 'State_Odd.dat'): # Only *.geom file.
        State_Odd_matches.append(path.join(file))
    if(len(State_Odd_matches) >= 2):
        print "WSC_Sim: there is more than one State_Odd.dat input file. there must be only one!"
    elif(len(State_Odd_matches) == 0):
        t_Odd = -1
    else: # len(State_Odd_matches) = 1
        State_Odd = open(State_Odd_matches[0],'r')
        State_Odd.readline().rstrip() # Title.
        State_Odd.readline().rstrip() # Starting time.
        ligne = State_Odd.readline().split() # t = time.
        t_Odd = float(ligne[0])

    return t_Even,t_Odd

def WSC_Sim(cwd):

    # cwd : Current working directory where this subroutine is called.

    # This subroutine runs several WSC simulations in using the input state files in order to avoid the memory leak problem.

    # Searching the *.in input file.
    INmatches = []
    files = [f for f in listdir('.') if path.isfile(f)] # Loop over the files of the cwd.
    for file in filter(files, '*.in'): # Only *.in file.
        INmatches.append(path.join(file))
    if(len(INmatches) >= 2):
        print "WSC_Sim: there is more than one *.in input file. there must be only one!"
        exit()

    # Searching the *.geom input file.
    GEOMmatches = []
    files = [f for f in listdir('.') if path.isfile(f)] # Loop over the files of the cwd.
    for file in filter(files, '*.geom'): # Only *.geom file.
        GEOMmatches.append(path.join(file))
    if(len(GEOMmatches) >= 2):
        print "WSC_Sim: there is more than one *.geom input file. there must be only one!"
        exit()

    # Reading the starting time, the time step and the number of time steps.
    InFile = open(INmatches[0],'r')
    for j in range(0,45): # Head lines.
        InFile.readline().rstrip()
    ligne = InFile.readline().split() # t0 = starting time.
    t0 = float(ligne[5])
    ligne = InFile.readline().split() # dt = time step.
    dt = float(ligne[5])
    ligne = InFile.readline().split() # Nt = Number of time steps.
    Nt = int(ligne[7])

    # Final time
    Tf = Nt*dt

    # WSC simulations until the final time is reached.
    niter = 0
    t = t0
    while (t < Tf):

        niter = niter + 1

        if(niter == 1):

            # Try to find a state input file.
            t_Even,t_Odd = Read_State(cwd)

            # Run the WSC code with or without a state input file.
            if(t_Even < 0 and t_Odd < 0): # No state input file at the beginning of the simulation (normal use of the WSC code).
                Job = Popen(["weakscatterer.exe",INmatches[0],GEOMmatches[0]]).wait()
            else: # There is already a state input file (maybe because the user had to stop the WSC code for any reason).
                if(t_Even > t_Odd): # State_Even.dat is used.
                    t = t_Even
                    # Saving the new input state file.
                    STATE_file = "State_"+str(t)+".dat"
                    copyfile("State_Even.dat", STATE_file)
                else: # State_Odd.dat is used.
                    t = t_Odd
                    # Saving the new input state file.
                    STATE_file = "State_"+str(t)+".dat"
                    copyfile("State_Odd.dat", STATE_file)
                Job = Popen(["weakscatterer.exe",INmatches[0],GEOMmatches[0],STATE_file]).wait()

        else: # niter != 1
            Job = Popen(["weakscatterer.exe",INmatches[0],GEOMmatches[0],STATE_file]).wait()

        # Reading the state input files.
        t_Even,t_Odd = Read_State(cwd)
        if(t_Even < 0 and t_Odd < 0):
            print "WSC_Sim: neither State_Even.dat not State_Odd.dat have been created."
            exit()
        if(t_Even < Nt and t_Odd < Nt): # If the simulation is not finished.
            if(t_Even > t_Odd): # State_Even.dat is used.
                t = t_Even
                # Saving the new input state file.
                STATE_file = "State_"+str(t)+".dat"
                copyfile("State_Even.dat", STATE_file)
            else: # State_Odd.dat is used.
                t = t_Odd
                # Saving the new input state file.
                STATE_file = "State_"+str(t)+".dat"
                copyfile("State_Odd.dat", STATE_file)

