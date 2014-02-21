#!/usr/bin/env python
#
# Raman off-resonant activity calculator
# using CRYSTAL as a back-end.
#
# Contributors: Alexandr Fonari (Georgia Tech)
# MIT license, 2013
#
#
def parse_fort34_header(fort34_fh):
    import sys
    from math import sqrt
    #
    fort34_fh.seek(0) # just in case
    #
    header = [fort34_fh.next() for x in range(5)] # read first 5 lines
    #
    vol = 0.0
    b = []
    #
    for i in range(1,4): #cartesian components of direct lattice vectors
        b.append( [float(s) for s in header[i].split()] )
    #
    vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
          b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]
    #
    symm_ops = int(header[-1])
    print "[parse_fort34_header]: Number of symmetry operations: %d" % symm_ops
    header.extend( [fort34_fh.next() for x in range(symm_ops*4)] )# each symm_op has 4 lines
    #
    header.extend( [fort34_fh.next()] ) # nat
    #
    nat_asym, nat_tot = [int(x) for x in header[-1].split()] # number of the irreducible atoms and total number of atoms in the primitive cell
    header[-1] = " %d\n" % nat_asym                          # I know, dirty hack
    #
    return nat_tot, vol, header
#
def parse_env_params(params):
    import sys
    #
    tmp = params.strip().split('_')
    if len(tmp) != 4:
        print "[parse_env_params]: ERROR there should be exactly four parameters"
        sys.exit(1)
    #
    [first, last, nderiv, step_size] = [int(tmp[0]), int(tmp[1]), int(tmp[2]), float(tmp[3])]
    #
    return first, last, nderiv, step_size
#
def get_modes_from_OUTCAR(outcar_fh, nat):
    import sys
    import re
    from math import sqrt
    eigvals    =  [ 0.0 for i in range(nat*3) ]
    eigvecs    =  [ 0.0 for i in range(nat*3) ]
    norms      =  [ 0.0 for i in range(nat*3) ]
    activity   =  [ '' for i in range(nat*3) ]
    atom_number = [ 0 for i in range(nat) ]
    pos        =  [ 0.0 for i in range(nat) ]
    asym        = [ 0.0 for i in range(nat) ]
    #
    outcar_fh.seek(0) # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "Eigenvectors after division by SQRT(mass)" in line:
            outcar_fh.readline() # empty line
            outcar_fh.readline() # Eigenvectors and eigenvalues of the dynamical matrix
            outcar_fh.readline() # ----------------------------------------------------
            outcar_fh.readline() # empty line
            #
            for i in range(nat*3): # all frequencies should be supplied, regardless of those requested to calculate
                outcar_fh.readline() # empty line
                p = re.search(r'^\s*(\d+).+?([-\.\d]+) cm-1 (\w)', outcar_fh.readline())
                eigvals[i] = float(p.group(2))
                activity[i] = p.group(3)
                #
                outcar_fh.readline() # X         Y         Z           dx          dy          dz
                eigvec = []
                #
                for j in range(nat):
                    tmp = outcar_fh.readline().split()
                    #
                    if i == 0: # get atomic positions only once
                        atom_number[j] = int(tmp[0])
                        pos[j]         = [ float(tmp[x]) for x in range(1,4) ]
                        asym[j]        = int(tmp[-1]) # is this atom in the asymmetric unit
                    #
                    eigvec.append([ float(tmp[x]) for x in range(4,7) ])
                    #
                eigvecs[i] = eigvec
                norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
            #
            return pos, asym, atom_number, eigvals, activity, eigvecs, norms
        #
    print "[get_modes_from_OUTCAR]: ERROR Couldn't find 'Eigenvectors after division by SQRT(mass)' in OUTCAR, exiting..."
    sys.exit(1)
#
###########################################
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False
###########################################
def get_epsilon_from_OUTCAR(outcar_fh):
    #
    eps = [[0.0 for i in range(3)] for j in range(3)]
    #
    outcar_fh.seek(0) # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "COMPONENT    ALPHA        EPSILON       CHI(1)" in line: # geeeeeez
            while True:
                line = outcar_fh.readline().strip()
                if not line: # empty line
                    break
                #
                direction, alpha, eps1, chi = line.split()
                for case in switch(direction.strip()):
                    if case('XX'):
                        eps[0][0] = float(eps1)
                        break
                    if case('XY'):
                        eps[0][1]= float(eps1)
                        break
                    if case('XZ'):
                        eps[0][2]= float(eps1)
                        break
                    if case('YY'):
                        eps[1][1]= float(eps1)
                        break
                    if case('YZ'):
                        eps[1][2]= float(eps1)
                        break
                    if case('ZZ'):
                        eps[2][2]= float(eps1)
                        break
                #
                #eps[2][0] = eps[0][2]
                #eps[1][0] = eps[0][1]
                #eps[2][1] = eps[1][2]
            #
            return eps
            break # while True
    #
    # no eps - no next mode
    raise RuntimeError("[get_epsilon_from_OUTCAR]: ERROR Couldn't find dielectric tensor in OUTCAR")
    return 1
#
if __name__ == '__main__':
    import sys
    from math import pi
    from shutil import move
    import os
    import datetime
    import time
    #import argparse
    import optparse
    #
    print ""
    print "    Raman off-resonant activity calculator,"
    print "    using CRYSTAL as a back-end."
    print ""
    print "    Contributors: Alexandr Fonari (Georgia Tech)"
    print "    MIT License, 2013"
    print "    URL: http://..."
    print "    Started at: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    print ""
    #
    description  = "Set environment variables:\n"
    description += "    CRYSTAL_RAMAN_RUN='runmpi09 INCAR'\n"
    description += "    CRYSTAL_RAMAN_PARAMS='[first-mode]_[last-mode]_[nderiv]_[step-size]'\n\n"
#    description += "One-liner bash is:\n"
#    description += "CRYSTAL_RAMAN_RUN='mpirun vasp' CRYSTAL_RAMAN_PARAMS='1 2 2 0.01' python vasp_raman.py -h"

    CRYSTAL_RAMAN_RUN = os.environ.get('CRYSTAL_RAMAN_RUN')
    if CRYSTAL_RAMAN_RUN == None:
        print "[__main__]: ERROR Set environment variable 'CRYSTAL_RAMAN_RUN'"
        print ""
        parser.print_help()
        sys.exit(1)
    print "[__main__]: CRYSTAL_RAMAN_RUN='"+CRYSTAL_RAMAN_RUN+"'"
    #
    CRYSTAL_RAMAN_PARAMS = os.environ.get('CRYSTAL_RAMAN_PARAMS')
    if CRYSTAL_RAMAN_PARAMS == None:
        print "[__main__]: ERROR Set environment variable 'CRYSTAL_RAMAN_PARAMS'"
        print ""
        parser.print_help()
        sys.exit(1)
    print "[__main__]: CRYSTAL_RAMAN_PARAMS='"+CRYSTAL_RAMAN_PARAMS+"'"
    #
    first, last, nderiv, step_size = parse_env_params(CRYSTAL_RAMAN_PARAMS)
    assert first >= 1,    '[__main__]: First mode should be equal or larger than 1'
    assert last >= first, '[__main__]: Last mode should be equal or larger than first mode'
    assert nderiv == 2,   '[__main__]: At this time, nderiv = 2 is the only supported'
    disps = [-1, 1]      # hardcoded for
    coeffs = [-0.5, 0.5] # three point stencil (nderiv=2)
    #
    try:
        fort34_fh = open('FORT34.phon', 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open input file FORT34.phon, exiting...\n"
        sys.exit(1)
    #
    nat, vol, fort34_header = parse_fort34_header(fort34_fh)
    fort34_fh.close()
    #
    try:
        outcar_fh = open('OUTCAR.phon', 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open OUTCAR.phon, exiting...\n"
        sys.exit(1)
    #
    pos, asym, atom_number, eigvals, activity, eigvecs, norms = get_modes_from_OUTCAR(outcar_fh, nat)
    outcar_fh.close()
    #
    output_fh = open('crystal_raman.dat', 'w')
    output_fh.write("# mode    freq(cm-1)    alpha    beta2    activity\n")
    for i in range(first-1, last):
        eigval = eigvals[i]
        eigvec = eigvecs[i]
        norm = norms[i]
        #
        print ""
        print "[__main__]: Mode #%i: frequency %10.7f cm-1; norm: %10.7f" % ( i+1, eigval, norm )
        #
        if activity[i] != 'A':
            print "[__main__]: Mode inactive, skipping..."
            continue
        #
        ra = [[0.0 for x in range(3)] for y in range(3)]
        for j in range(len(disps)):
            disp_filename = 'OUTCAR.%04d.%+d.out' % (i+1, disps[j])
            #
            try:
                outcar_fh = open(disp_filename, 'r')
                print "[__main__]: File "+disp_filename+" exists, parsing..."
            except IOError:
                print "[__main__]: File "+disp_filename+" not found, preparing displaced INCAR.gui"
                poscar_fh = open('INCAR.gui', 'w')
                poscar_fh.write("".join(fort34_header))
                #
                for k in range(nat): # do the deed
                    if asym[k] == 0: continue # this atom is NOT in the asymmetric unit, skip!
                    #
                    pos_disp = [ pos[k][l] + eigvec[k][l]*step_size*disps[j]/norm for l in range(3)]
                    poscar_fh.write( "%3d %15.10f %15.10f %15.10f\n" % (atom_number[k], pos_disp[0], pos_disp[1], pos_disp[2]) )
                    #print '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f' % (pos[k][0], pos[k][1], pos[k][2], dis[k][0], dis[k][1], dis[k][2])
                poscar_fh.close()
                #
                # run CRYSTAL here
                print "[__main__]: Running CRYSTAL..."
                os.system(CRYSTAL_RAMAN_RUN)
                try:
                    move('INCAR.out', disp_filename)
                except IOError:
                    print "[__main__]: ERROR Couldn't find INCAR.out file, exiting..."
                    sys.exit(1)
                #
                outcar_fh = open(disp_filename, 'r')
            #
            try:
                eps = get_epsilon_from_OUTCAR(outcar_fh)
                outcar_fh.close()
            except Exception, err:
                print err
                print "[__main__]: Moving "+disp_filename+" back to 'OUTCAR' and exiting..."
                move(disp_filename, 'OUTCAR')
                sys.exit(1)
            #
            for m in range(3):
                for n in range(3):
                    ra[m][n]   += eps[m][n] * coeffs[j]/step_size * norm * vol/(4.0*pi)
            #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3
        #
        alpha = (ra[0][0] + ra[1][1] + ra[2][2])/3.0
        beta2 = ( (ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 + (ra[1][1] - ra[2][2])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2) )/2.0
        print ""
        print "! %4i  freq: %10.5f  alpha: %10.7f  beta2: %10.7f  activity: %10.7f " % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2)
        output_fh.write("%i  %10.5f  %10.7f  %10.7f  %10.7f\n" % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2))
        output_fh.flush()
    #
    output_fh.close()
    sys.exit(0)
    # done.
