#!/usr/bin/env python

from __future__ import division
import os
import re
import sys
import numpy as np
import argparse as arg
from elements import ELEMENTS

#au2ang = 0.5291771


def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''Generates .xyz trajectories for
                                each Normal Mode from a QM frequency calculation''',
                                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input files
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('-f', '--filename',
                     default=None, type=str, dest="File", required=True,
                     help='''Excited State file''')

    inp.add_argument('-s', '--sel', default=None, nargs='+', type=str,
                     dest='AtomSel', help='''Atom Selection.''')

    inp.add_argument('--scale', default=5, type=int,
                     dest='Scale', help='''Scaling factor for displacement
                     vectors in the VMD script.''')

    #
    # Output files
    #
    out = parser.add_argument_group("Output Data")

    out.add_argument('-o', '--outdir', default="normal_modes",
                     type=str, dest="OutDir", help='''Output folder''')

    #
    # Parse and create the Options Dictionary
    #
    args = parser.parse_args()
    Opts = vars(args)

    if Opts['AtomSel']:
        Opts['AtomSel'] = read_sel(Opts['AtomSel'])

    return Opts


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  list(map(int, to_extend))
            extended += range(sel[0],sel[1]+1,1)

        else:

            extended.append(int(idx))

    return extended


def read_sel(string):

    string =  ','.join(string).replace(',,',',')

    try:
        f = open(string, 'r')
        string = f.readlines()
        f.close()
        string =  ','.join(string).replace(',,',',')
        string = string.replace(',', ' ')
        string = list(map(lambda x: x - 1, extend_compact_list(string)))

    except IOError:
        string = string.replace(',', ' ')
        string = list(map(lambda x: x - 1, extend_compact_list(string)))

    return string


def parsefreqs_QChem(filename):

    with open(filename) as f:

        for line in f:

            #
            # Geometry
            #
            if "Standard Nuclear Orientation" in line:

                # This next line guarantees we only retrieve the last structure
                # in case a double job opt+freq is contained in the logfile
                structure = []
                line = skiplines(f, 2)
                data = line.split()

                while len(data) == 5:

                    atom = data[1]
                    atom_x = float(data[2])
                    atom_y = float(data[3])
                    atom_z = float(data[4])
                    structure.append([atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

            #
            # Vibrational Properties
            #
            if "VIBRATIONAL ANALYSIS" in line:

                freqs = []
                forcecns = []
                redmasses = []
                modes = []
                line = skiplines(f, 10)

                while line:

                    #
                    # Frequencies
                    #
                    if "Frequency" in line:
                        tmpfreqs = list(map(float, line.split()[1:]))
                        freqs.extend(tmpfreqs)

                    #
                    # Force Constants
                    #
                    if "Force Cnst" in line:
                        tmpcns = list(map(float, line.split()[2:]))
                        forcecns.extend(tmpcns)

                    #
                    # Reduced Masses
                    #
                    if "Red. Mass" in line:
                        tmpredms = list(map(float, line.split()[2:]))
                        redmasses.extend(tmpredms)

                    #
                    # Cartesian Displacements
                    #
                    if "Raman Active" in line:
                        line = skiplines(f, 1)
                        disps = []

                        while "TransDip" not in line:
                            data = list(map(float, line.split()[1:]))
                            N = len(data) // 3

                            if not disps:
                                for n in range(N):
                                    disps.append([])

                            for n in range(N):
                                disps[n].append(data[3*n:(3*n)+3])

                            line = skiplines(f)

                        modes.extend(disps)

                    try:
                        line = skiplines(f)

                    except StopIteration:
                        break


        Z_atoms = [ ELEMENTS[x[0]].number for x in structure ]
        masses = np.array([ ELEMENTS[x[0]].mass for x in structure ])
        atoms = [ x[0] for x in structure ]
        coords = np.array([ x[1:] for x in structure ])
        freqs = np.array(freqs)
        forcecns = np.array(forcecns)
        redmasses = np.array(redmasses)
        modes = np.array(modes)

        return Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, modes


def parsefreqs_G09(filename):

    with open(filename) as f:

        structures = []
        freqs = []
        forcecns = []
        redmasses = []
        masses = []
        modes = []
        HPmodes = False

        struct_done = False

        for line in f:

            #
            # Atomic Masses
            #
            if "AtmWgt" in line and not struct_done:
                data = list(map(float, line.split()[1:]))
                masses.extend(data)

            #
            # Non reoriented geometry
            #
            if "Input orientation" in line:
                structure = []
                line = skiplines(f, 4)
                data = line.split()

                while len(data) == 6:

                    Z_atom = int(data[1])
                    atom_x = float(data[3])
                    atom_y = float(data[4])
                    atom_z = float(data[5])
                    structure.append([Z_atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                structures.append(structure)
                struct_done = True

            #
            # Reoriented geometry, overwrite previous
            #
            if "Standard orientation" in line:
                structure = []
                line = skiplines(f, 4)
                data = line.split()

                while len(data) == 6:

                    Z_atom = int(data[1])
                    atom_x = float(data[3])
                    atom_y = float(data[4])
                    atom_z = float(data[5])
                    structure.append([Z_atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                structures.append(structure)
                struct_done = True

            #
            # Vibrational Properties
            #

            #
            # HPmodes
            #
            if "Frequencies ---" in line:
                HPmodes = True
                NAtoms = len(structure)
                tmpfreqs = list(map(float, line.split()[2:]))
                freqs.extend(tmpfreqs)

            if "Reduced masses ---" in line:
                tmpredms = list(map(float, line.split()[3:]))
                redmasses.extend(tmpredms)

            if "Force constants ---" in line:
                tmpcns = list(map(float, line.split()[3:]))
                forcecns.extend(tmpcns)

            if "Coord Atom Element" in line:
                disps = []

                for n in range(3 * NAtoms):
                    line = skiplines(f)
                    data = line.split()
                    atomindex = int(data[1]) - 1
                    numbers = list(map(float, data[3:]))
                    numbermodes = len(numbers)

                    if not disps:
                        for mode in range(numbermodes):
                            disps.append([[] for x in range(0, NAtoms)])

                    for mode in range(numbermodes):
                        disps[mode][atomindex].append(numbers[mode])

                modes.extend(disps)

            #
            # No HPmodes
            #
            if not HPmodes:

                if "Frequencies --" in line:
                    NAtoms = len(structure)
                    tmpfreqs = list(map(float, line.split()[2:]))
                    freqs.extend(tmpfreqs)

                if "Red. masses --" in line:
                    tmpredms = list(map(float, line.split()[3:]))
                    redmasses.extend(tmpredms)

                if "Frc consts --" in line:
                    tmpcns = list(map(float, line.split()[3:]))
                    forcecns.extend(tmpcns)

                if "X      Y      Z" in line:
                    line = skiplines(f)
                    disps = []

                    i = 1
                    while i <= NAtoms:
                        data = list(map(float, line.split()[2:]))
                        N = len(data) // 3

                        if not disps:
                            for n in range(N):
                                disps.append([])

                        for n in range(N):
                            disps[n].append(data[3*n:(3*n)+3])

                        line = skiplines(f)
                        i += 1

                    modes.extend(disps)

        if not masses:
            masses = [ ELEMENTS[x[0]].mass for x in structure ]

        Z_atoms = [ int(x[0]) for x in structure ]
        masses = np.array(masses)
        atoms = [ ELEMENTS[x[0]].symbol for x in structure ]
        coords = np.array([ x[1:] for x in structures[-1] ])
        freqs = np.array(freqs)
        forcecns = np.array(forcecns)
        redmasses = np.array(redmasses)
        modes = np.array(modes)

        return Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, modes


def guess(filename):
    '''Returns the correct class needed to parse filename, if it exists.'''

    #
    # Dictionary of unique sentences in QM packages output files to guess
    # the correct parser to use
    #
    filetypes = {}

    filetypes[" This is part of the Gaussian(R) 16 program."] = "G09"
    filetypes["A Quantum Leap Into The Future Of Chemistry"] = "QChem"

    filetype = None
    done = False
    with open(filename) as f:

        for line in f:
            for sentence in filetypes.keys():

                if sentence in line:
                    filetype = filetypes[sentence]
                    done = True
                    break

            # once the type has been identified, exit the cycles
            if done:
                break

    if not filetype:
        print(" %s" % filename)
        print(" File type not known")
        sys.exit()

    return filetype


def save_visdisps(coords, disps, sel=None, filename="structure"):

    header = ("#\n"
              "# VMD script to draw vectors\n"
              "#\n"
              "menu main on\n"
              "display projection orthographic\n"
              "display depthcue off\n"
              "display nearclip set 0.01\n"
              "axes location lowerleft\n"
              "\n"
              "#\n"
              "# VMD functions to draw a vector\n"
              "#\n"
              "proc vmd_draw_arrow {mol start end} {\n"
              "    set length [veclength [vecsub $end $start]]\n"
              "    set conelen [expr max(0.4,0.2*$length) ]\n"
              "    set scale [expr max(0.5,(1.0-$conelen/$length))]\n"
              "\n"
              "    set middle [vecadd $start [vecscale $scale [vecsub $end $start]]]\n"
              "    graphics $mol cylinder $start $middle radius 0.05\n"
              "    puts [list cone $middle $end radius 0.15]\n"
              "    graphics $mol cone $middle $end radius 0.15\n"
              "}\n"
              "\n"
              "proc vmd_draw_vector { mol pos val } {\n"
              "    set end   [ vecadd $pos [ vecscale +1 $val ] ]\n"
              "    vmd_draw_arrow $mol $pos $end\n"
              "}\n"
              "\n"
              "\n")

    #
    # Save the geometry and write the VMD script file with transition dips
    #
    with open('%s.vmd' % filename, 'w') as f:
        f.write(header)
        f.write("mol new equilibrium.xyz type xyz\n")
        f.write("mol new %s.xyz type xyz\n" % os.path.split(filename)[-1])

        if sel:
            f.write("mol rep Licorice 0.2 10\n")
            f.write("mol selection index %d to %d\n" % (sel[0], sel[-1]))
            f.write("mol addrep top\n")

        f.write("mol showrep 0 0 off\n")
        f.write("\n")

        for N, coord in enumerate(coords[sel][0]):

            disp = disps[N] * Opts['Scale']
            if np.linalg.norm(disp) > 0.0015 * Opts['Scale']:
                cmd = "graphics 0 color green; vmd_draw_vector 0 {%8.4f %8.4f %8.4f} {%8.4f %8.4f %8.4f}\n"
                data = [coord[0], coord[1], coord[2], disp[0], disp[1], disp[2]]
                f.write(cmd % tuple(data))

            # scale = Opts['Scale'] / np.linalg.norm(disps[N])
            # disp = disps[N] * scale
            # cmd = "graphics 0 color green; vmd_draw_vector 0 {%8.4f %8.4f %8.4f} {%8.4f %8.4f %8.4f}\n"
            # data = [coord[0], coord[1], coord[2], disp[0], disp[1], disp[2]]
            # f.write(cmd % tuple(data))

    return


if __name__ == '__main__':

    Opts = options()

    #
    # Choose the parser depending on the method and on the input files
    #
    parsertype = guess(Opts["File"])

    if parsertype == "G09":
        freqparser = parsefreqs_G09

    elif parsertype == "QChem":
        freqparser = parsefreqs_QChem

    Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, xdisps = freqparser(Opts['File'])

    #
    # Create OutDir
    #
    if not os.path.exists(Opts['OutDir']):
        os.makedirs(Opts['OutDir'])

    #
    # Save Equilibrium Geometry
    #
    eqfile = os.path.join(Opts['OutDir'], "equilibrium")
    with open("%s.xyz" % eqfile, "w") as f:
        f.write("%d\n\n" % len(Z_atoms))
        np.savetxt(f, np.c_[Z_atoms, coords], fmt="%3d %12.6f %12.6f %12.6f")

    #
    # Cycle over Normal Modes
    #
    for N, xdisp in enumerate(xdisps, start=1):

        mode_anim = []

        #
        # Displace the eq geometry
        #
        for i in np.linspace(-1, 1, 20):

            coor_disp = coords + xdisp * i
            disp_geom = np.c_[Z_atoms, coor_disp]
            mode_anim.append(disp_geom)

        mode_anim = np.array(mode_anim)
        outfile = os.path.join(Opts['OutDir'], "mode_%03d_%d" % (N, freqs[N-1]))

        #
        # Save movie
        #
        with open("%s.xyz" % outfile, "w") as f:

            for step in mode_anim:
                f.write("%d\n\n" % len(Z_atoms))
                np.savetxt(f, step, fmt="%3d %12.6f %12.6f %12.6f")

            for step in mode_anim[::-1]:
                f.write("%d\n\n" % len(Z_atoms))
                np.savetxt(f, step, fmt="%3d %12.6f %12.6f %12.6f")

        save_visdisps(coords, xdisp, sel=Opts['AtomSel'], filename="%s" % outfile)
