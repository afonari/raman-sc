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
    N = 10 # header lines number
    header = [fort34_fh.next() for x in range(N)]
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
    num_symm = int(header[4]) # number of symmetry operators
    if num_symm != 1:
        print "[parse_fort34]: Number of symmetry operations should be 1 (line 5 in GUI file), exiting..."
        sys.exit(1)
    #
    nat = int(header[9]) # number of irreducible atoms in the primitive cell
    #
    return nat, vol, header
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
    eigvals = [ 0.0 for i in range(nat*3) ]
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms   = [ 0.0 for i in range(nat*3) ]
    pos     = [ 0.0 for i in range(nat) ]
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
                p = re.search(r'^\s*(\d+).+?([\.\d]+) cm-1', outcar_fh.readline())
                eigvals[i] = float(p.group(2))
                #
                outcar_fh.readline() # X         Y         Z           dx          dy          dz
                eigvec = []
                #
                for j in range(nat):
                    tmp = outcar_fh.readline().split()
                    if tmp[-1] == 0: continue # in asymmetric unit or not
                    if i == 0: pos[j] = [ float(tmp[x]) for x in range(3) ] # get atomic positions only once
                    #
                    eigvec.append([ float(tmp[x]) for x in range(3,6) ])
                    #
                eigvecs[i] = eigvec
                norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
            #
            return pos, eigvals, eigvecs, norms
        #
    print "[get_modes_from_OUTCAR]: ERROR Couldn't find 'Eigenvectors after division by SQRT(mass)' in OUTCAR. Use 'NWRITE=3' in INCAR. Exiting..."
    sys.exit(1)
#
def parse_outcar(outcar_fn):
    import re
    import sys
    eps = [[0.0 for i in range(3)] for j in range(3)]

    outcar_fh = open(outcar_fn, 'r')
    while True:
        line = outcar_fh.readline()
        if not line:
            break

        if "COMPONENT    ALPHA        EPSILON       CHI(1)" in line:
            
                eps[0][0] = float(outcar_fh.readline()[21:42])
                eps[0][1] = float(outcar_fh.readline()[21:42])
                eps[0][2] = float(outcar_fh.readline()[21:42])
                eps[1][1] = float(outcar_fh.readline()[21:42])
                eps[1][2] = float(outcar_fh.readline()[21:42])
                eps[2][2] = float(outcar_fh.readline()[21:42])
                eps[2][0] = eps[0][2]
                eps[1][0] = eps[0][1]
                eps[2][1] = eps[1][2]
                return eps
    #
    # no eps - no next mode
    print "Couldn't find dielectric tensor in "+outcar_fn+", exiting..."
    sys.exit(1)

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
    pos, eigvals, eigvecs, norms = get_modes_from_OUTCAR(outcar_fh, nat)
    outcar_fh.close()

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
        ra = [[0.0 for x in range(3)] for y in range(3)]
        for j in range(len(disps)):
            disp_filename = 'OUTCAR.%04d.%+d.out' % (i+1, disps[j])
            #
            try:
                outcar_fh = open(disp_filename, 'r')
                print "[__main__]: File "+disp_filename+" exists, parsing..."
            except IOError:
                print "[__main__]: File "+disp_filename+" not found, preparing displaced POSCAR"
                poscar_fh = open('fort34.in', 'w')
                poscar_fh.write("\n".join(fort34_header))
                #
                for k in range(nat):
                    pos_disp = [ pos[k][l] + eigvec[k][l]*step_size*disps[j]/norm for l in range(3)]
                    poscar_fh.write( "%15.10f %15.10f %15.10f\n" % (pos_disp[0], pos_disp[1], pos_disp[2]) )
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
#    #
#    # 2nd: open fort34 file
#    try:
#        fort34_fn = sys.argv[1]
#        fort34_fh = open(fort34_fn, 'r')
#    except IndexError:
#        print "Run code as:\n"
#        print "    CRYSTAL-Raman.py fort34.in\n"
#        sys.exit("No argument found, exiting...\n")
#    except IOError:
#        sys.exit("Couldn't open "+fort34_fn+" file, exiting...\n")
#    #
#    # 3rd: everything else
#    print "Parsing "+fort34_fn+"...\n"
#    first, last, nderiv, step_size, input_fn, modes_fn, am_fn, nat, vol, pos, atomic_numbers, fort34_header = parse_fort34(fort34_fh)
#    #
#    eig_fh = open(modes_fn, 'r')
#    #
#    if nderiv == 2:
#        disps = [-1, 1]
#        coeffs = [-0.5, 0.5]
#    else:
#        print "Unknown value for NDERIV (use 2), exiting..."
#        sys.exit(0)
#    #
#    try:
#        active_modes = [line.strip() for line in open(am_fn, 'r')]
#    except IOError:
#        sys.exit("Couldn't open "+am_fn+" file, exiting...\n")
#    #
#    print "Number of ions: %i, volume: %7.5f A^3" % (nat, vol)
#    print "Modes to be computed: %i - %i" % (first, last)
#    print "Derivation number and step size (A): %i, %7.5f" % (nderiv, step_size)
#    print "Input filename: %s" % input_fn
#    print "Modes filename: %s" % modes_fn
#    print "Active modes filename: %s" % am_fn
#    print "CRYSTAL invoke command: %s" % CRYSTAL_RUN 
#    print ""
#    #
#    OUT_FN=input_fn+'.out'
#    GUI_FN=input_fn+'.gui'
#    #
#    eigvals, eigvecs, norms = read_qe_dynmat(eig_fh, nat)
#    #
#    for i in range(first-1, last):
#        eigval = eigvals[i]
#        eigvec = eigvecs[i]
#        norm = norms[i]
#    
#        print "Mode #: %i, frequency: %10.7f cm-1, norm: %10.7f" % ( i+1, eigval, norm )
#        if active_modes[i] == 'I':
#            print "Mode #: %i - raman inactive, skipping..." % (i+1)
#            continue
#
#        ra = [[0.0 for x in range(3)] for y in range(3)]
#        for j in range(len(disps)):
#            disp_filename = 'OUTCAR.%04d.%+d.out' % (i+1, disps[j])
#            #
#            try:
#                with open(disp_filename, 'r'):
#                    file_exists =True
#                    pass
#            except IOError:
#                file_exists=False
#            #
#            if file_exists:
#                print "File "+disp_filename+" exists, parsing..."
#            else:
#                print "File "+disp_filename+" not found, running CRYSTAL...\n"
#                gui_fh = open(GUI_FN, 'w')
#                gui_fh.write(fort34_header)
#                #
#                for k in range(nat):
#                    pos_disp = [ pos[k][l] + eigvec[k][l]*step_size*disps[j]/norm for l in range(3)]
#                    gui_fh.write( '%5d     %15.12f     %15.12f     %15.12f\n' % (atomic_numbers[k], pos_disp[0], pos_disp[1], pos_disp[2]) )
#                    #print '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f' % (pos[k][0], pos[k][1], pos[k][2], dis[k][0], dis[k][1], dis[k][2])
#                gui_fh.close()
#                #
#                #run CRYSTAL here
#                os.system(CRYSTAL_RUN)
#                #
#                try:
#                    move(OUT_FN, disp_filename)
#                except IOError:
#                    sys.exit("Couldn't open "+OUT_FN+" file, exiting...\n")
#            #
#            print "Parsing "+disp_filename+"...\n"
#            eps = parse_outcar(disp_filename)
#            print "Found eps in "+disp_filename+":"
#            print eps
#            #
#            for m in range(3):
#                for n in range(3):
#                    ra[m][n]   += eps[m][n] * coeffs[j]/step_size * norm * vol/(4.0*pi)
#            #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3
#        #
#        alpha = (ra[0][0] + ra[1][1] + ra[2][2])/3.0
#        beta2 = ( (ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 + (ra[1][1] - ra[2][2])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2) )/2.0
#        print 'for mode %i: %10.5f cm-1; alpha: %10.7f ; beta2: %10.7f ; activity: %10.7f ' % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2)
