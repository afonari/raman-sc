#!/usr/bin/env python
#
# Helper script to generate Eigenvectors/Eigenvalues from CRYSTAL output
#
# Contributors: Alexandr Fonari (Georgia Tech)
# MIT license, 2013
#
#
from __future__ import division
#
def jacobi(ainput):
    # from NWChem/contrib/python/mathutil.py
    # possible need to rewrite due to licensing issues
    #
    from math import sqrt
    #
    a = [[ ainput[i][j] for i in range(len( ainput[j] )) ] for j in range(len( ainput )) ] # copymatrix
    n = len(a)
    m = len(a[0])
    if n != m:
        raise 'jacobi: Matrix must be square'
    #
    for i in range(n):
        for j in range(m):
            if a[i][j] != a[j][i]:
                raise 'jacobi: Matrix must be symmetric'
    #
    tolmin = 1e-14
    tol = 1e-4
    #
    v = [[0.0 for i in range(n)] for j in range(n)] # zeromatrix
    for i in range(n):
        v[i][i] = 1.0
    #
    maxd = 0.0
    for i in range(n):
        maxd = max(abs(a[i][i]),maxd)
    #
    for iter in range(50):
        nrot = 0
        for i in range(n):
            for j in range(i+1,n):
                aii = a[i][i]
                ajj = a[j][j]
                daij = abs(a[i][j])
                if daij > tol*maxd: # Screen small elements
                    nrot = nrot + 1
                    s = aii - ajj
                    ds = abs(s)
                    if daij > (tolmin*ds): # Check for sufficient precision
                        if (tol*daij) > ds:
                            c = s = 1/sqrt(2.)
                        else:
                            t = a[i][j]/s
                            u = 0.25/sqrt(0.25+t*t)
                            c = sqrt(0.5+u)
                            s = 2.*t*u/c
                        #
                        for k in range(n):
                            u = a[i][k]
                            t = a[j][k]
                            a[i][k] = s*t + c*u
                            a[j][k] = c*t - s*u
                        #
                        for k in range(n):
                            u = a[k][i]
                            t = a[k][j]
                            a[k][i] = s*t + c*u
                            a[k][j]= c*t - s*u
                        #
                        for k in range(n):
                            u = v[i][k]
                            t = v[j][k]
                            v[i][k] = s*t + c*u
                            v[j][k] = c*t - s*u
                        #
                        a[j][i] = a[i][j] = 0.0
                        maxd = max(maxd,abs(a[i][i]),abs(a[j][j]))
        #
        if nrot == 0 and tol <= tolmin:
            break
        tol = max(tolmin,tol*0.99e-2)
    #
    if nrot != 0:
        print 'jacobi: [WARNING] Jacobi iteration did not converge in 50 passes!'
    #
    # Sort eigenvectors and values into increasing order
    e = [0.0 for i in range(n)] # zerovector
    for i in range(n):
        e[i] = a[i][i]
        for j in range(i):
            if e[j] > e[i]:
                (e[i],e[j]) = (e[j],e[i])
                (v[i],v[j]) = (v[j],v[i])
    #
    return (v,e)
#
def get_modes_from_HESSFREQ(hessfreq_fh, nat, masses):
    from math import sqrt
    #
    ELECTRONMASS_SI  = 9.10938215E-31  # Kg
    AMU_SI           = 1.660538782E-27 # Kg
    AMU_AU           = AMU_SI / ELECTRONMASS_SI
    AU_CMM1          = 0.2194746E+06
    dyn = []
    fc = [[0.0 for i in range(nat*3)] for j in range(nat*3)]
    freqs = []
    #
    hessfreq_fh.seek(0) # just in case
    #
    while True:
        line = hessfreq_fh.readline()
        if not line:
            break
        #
        tmp = line.strip().split()
        dyn.extend( [float(x) for x in line.strip().split()] )
    #
    for i in range(nat*3):
        for j in range(nat*3):
            t = i + j*nat*3
            fc[i][j] = dyn[i + j*nat*3]
    #
    for i in range(nat*3):
        for j in range(nat*3):
            fc[i][j] = fc[i][j]/(AMU_AU*sqrt(masses[i//3]*masses[j//3]))
    #
    for i in range(nat*3):
        for j in range(i):
            fc[j][i] = fc[i][j] # using the lower triangle... for some reasons ;)
    #
    eigvecs, eigvals = jacobi(fc)
#    import numpy as np
#    eigvals, eigvecs_ = np.linalg.eigh(np.array(fc))
#    for
    #
    for eigval in eigvals:
        freq = sqrt(abs(eigval))
        if(eigval < 0): freq = -freq
        #
        freqs.append(freq * AU_CMM1)
    #
    for i in range(nat*3):
        eigvecs[i] = [x/sqrt(masses[i//3]) for x in eigvecs[i]]
    #
    return freqs, eigvecs
#
def dump_OUTCAR_phon(outcar_fh, nat, pos, asym, freqs, eigvecs):
    outcar_fh.seek(0) # just in case
    outcar_fh.write("Eigenvectors after division by SQRT(mass)\n")
    outcar_fh.write("\n")
    outcar_fh.write("Eigenvectors and eigenvalues of the dynamical matrix\n")
    outcar_fh.write("----------------------------------------------------\n")
    outcar_fh.write("\n")
    #
    for i in range(nat*3):
        freq = freqs[i]
        evec = eigvecs[i]
        #
        outcar_fh.write("\n")
        outcar_fh.write("%5d  %5.6f cm-1\n" % (i+1, freq))
        outcar_fh.write(" X         Y         Z           dx          dy          dz\n")
        #
        for j in range(nat):
            a = int(asym[j])
            outcar_fh.write("%10.6f  %10.6f  %10.6f    %10.6f  %10.6f  %10.6f %d\n" % (pos[j][0], pos[j][1], pos[j][1], evec[3*j], evec[3*j+1], evec[3*j+2], a))

def parse_OUTCAR(outcar_fh):
    import sys
    nat = 0
    masses = []
    asym = []
    pos = []
    #
    outcar_fh.seek(0) # just in case
    #
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "ATOMS IN THE UNIT CELL:" in line: # ATOMS IN THE ASYMMETRIC UNIT   18 - ATOMS IN THE UNIT CELL:  140
            tmp = line.strip().split()
            nat = int(tmp[-1])
            #
            outcar_fh.readline() # ATOM              X/A                 Y/B                 Z/C
            outcar_fh.readline() # *******************************************************************************
            #
            for i in range(nat): # 1 T   6 C    -5.326841882804E-05 -1.351394382997E-04  1.224852639589E-01
                tmp = outcar_fh.readline().strip().split()
                #
                if tmp[1] == 'T':
                    asym.append(True)
                else:
                    asym.append(False)
        #
        if "CARTESIAN COORDINATES" in line: # CARTESIAN COORDINATES - PRIMITIVE CELL
            outcar_fh.readline() # *******************************************************************************
            outcar_fh.readline() # *      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)
            outcar_fh.readline() # *******************************************************************************
            #
            for i in range(nat):
                tmp = outcar_fh.readline().strip().split()
                pos.append( [ float(tmp[x]) for x in range(3,6) ] )
            #
        if "ATOMS ISOTOPIC MASS (AMU) FOR FREQUENCY CALCULATION" in line:
            outcar_fh.readline() # empty line
            while True:
                line = outcar_fh.readline()
                if not line: #EOF
                    break
                #
                if not line.strip(): # empty line
                    break
                #
                tmp = line.strip().split()
                for i in range(2,len(tmp),3): # fancy way to parse masses
                    masses.append( float(tmp[i]) )
    #
    assert nat == len(masses), "[parse_OUTCAR]: Number of atoms should be equal to the size of the 'masses' array"
    assert nat == len(asym),   "[parse_OUTCAR]: Number of atoms should be equal to the size of the 'asym' array"
    return nat, masses, asym, pos
#
if __name__ == '__main__':
    import sys
    #
    if len(sys.argv) < 3:
        print "[__main__]: ERROR Run as:"
        print "             crystal_modes_gen.py <freq.out> <HESSFREQ.dat>"
        print ""
        sys.exit(1)
    #
    try:
        outcar_fh = open(sys.argv[1], 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open first-argument file: '"+sys.argv[1]+"', exiting...\n"
        sys.exit(1)
    #
    nat, masses, asym, pos = parse_OUTCAR(outcar_fh)
    outcar_fh.close()
    #
    try:
        hessfreq_fh = open(sys.argv[2], 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open second-argument file: '"+sys.argv[2]+"', exiting...\n"
        sys.exit(1)
    #
    freqs, eigvecs = get_modes_from_HESSFREQ(hessfreq_fh, nat, masses)
    hessfreq_fh.close()
    #
    try:
        outcar_fh = open('OUTCAR.phon', 'w')
    except IOError:
        print "[__main__]: ERROR Couldn't open OUTCAR.phon for writing, exiting...\n"
        sys.exit(1)
    #
    dump_OUTCAR_phon(outcar_fh, nat, pos, asym, freqs, eigvecs)
    outcar_fh.close()
    print "'OUTCAR.phon' generated successfully!"
    print ""
    #
