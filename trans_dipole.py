import cclib
import sys
import numpy as np
import os
import re
import subprocess


# Functions needed to do post-processing analysis

def gauss_trans_dip_mom(name):
    """
    Returns a numpy array with the coordinates of the transition dipole moment 
    computed from a gaussian calculation.

    Parameters:
    -----------
    coord : : class 'numpy.ndarray'
    Numpy array with coordinates obtained from the
    log file of a gaussian calculation.

    mass : : class 'numpy.ndarray'
    Numpy array with masses obtained from the
    log file of a gaussian calculation.

    Returns: 
    --------
    :class: 'numpy.ndarray'
    Coordinates of the selected transition dipole moment
    """ 
    trs_dip_mom = []
    canprintlines = False
    with open(name) as file:
        for line in file:
          if line.startswith(' Ground to excited state transition electric dipole moments (Au):'):
              canprintlines = True
          elif line.startswith(' Ground to excited state transition velocity dipole moments'):
              canprintlines = False
          if canprintlines:
              trs_dip_mom.append(line)

    trs_dip_mom=np.array(trs_dip_mom)
    indi_dip_mom = []
    for i in range(1,len(trs_dip_mom)):
        indi_dip_mom.append(trs_dip_mom[i].split())

    indi_dip_mom = np.array(indi_dip_mom)
    test = []
    for i in range(len(indi_dip_mom)):
        test.append([indi_dip_mom[i][0], indi_dip_mom[i][1], indi_dip_mom[i][2],
                    indi_dip_mom[i][3], indi_dip_mom[i][4], indi_dip_mom[i][5]])

    test = np.array(test)
    dip_mom_all = test[1:,1:4]
    return dip_mom_all

def plot_dip_arrow(center_of_mass, dipole_moment,output):
   #####
   # Function to plot the dipole moments from gaussian
   # with the help of vmd and is done via vmd -e
   # output.vmd
   #####   
   with open('{0:s}.vmd'.format(output), 'w') as f:
        f.write("mol new {0:s}.xyz type xyz\n".format(output))
        f.write("draw material Opaque\n")
        f.write("draw color green\n")
        f.write("draw cylinder {%7.4f %7.4f %7.4f} {%7.4f %7.4f %7.4f} radius 0.1\n" %
                (center_of_mass[0], center_of_mass[1], center_of_mass[2], float(dipole_moment[0]), float(dipole_moment[1]), float(dipole_moment[2])))
   print ("{0:s}.vmd file has been created".format(output))

def com_coord(coord, masses):
  """
  Returns a numpy array with the coordinates of the center of 
  mass of a molecule computed from a gaussian calculation.

  Parameters:
  -----------
  coord : : class 'numpy.ndarray'
        Numpy array with coordinates obtained from the
        log file of a gaussian calculation.

  mass : : class 'numpy.ndarray'
        Numpy array with masses obtained from the
        log file of a gaussian calculation.
  Returns: 
  --------
  :class: 'numpy.ndarray'
      Coordinates of the center of mass of the computed molecule
  """ 
  x = []; y=[]; z=[]; mass=[];
  for i in range(natoms):
     x.append(coord[0][i][0])
     y.append(coord[0][i][1])
     z.append(coord[0][i][2])
     mass.append(masses[i])

  total_mass = np.sum(mass)
  x_com = [] ; y_com = [] ; z_com = [] ; com = []
  for i in range(natoms):                                 
    x_com.append(float(mass[i]*x[i]/total_mass))
    y_com.append(float(mass[i]*y[i]/total_mass))
    z_com.append(float(mass[i]*z[i]/total_mass))
    com = [np.sum(x_com), np.sum(y_com), np.sum(z_com)]
         
  com = np.array(com)
  return com

def get_xyz(ifile):
   #####
   # Function to plot the dipole moments from gaussian
   # with the help of vmd and is done via vmd -e
   # output.vmd
   #####   
   convert = subprocess.Popen(['ccwrite xyz %s' % ifile], shell=True)
   outs, errs = convert.communicate(timeout=120)
   print ("{}.xyz succesfully created".format(ifile))

def gauss_trans_mag_mom(name):
    """
    Returns a numpy array with the coordinates of the transition magnetic moment 
    computed from a gaussian calculation.

    Parameters:
    -----------
    coord : : class 'numpy.ndarray'
    Numpy array with coordinates obtained from the
    log file of a gaussian calculation.

    mass : : class 'numpy.ndarray'
    Numpy array with masses obtained from the
    log file of a gaussian calculation.

    Returns: 
    --------
    :class: 'numpy.ndarray'
      Coordinates of the selected magnetic dipole moment
    """ 
    trs_mag_mom = []
    canprintlines = False
    with open(name) as file:
        for line in file:
          if line.startswith(' Ground to excited state transition magnetic dipole moments (Au):'):
              canprintlines = True
          elif line.startswith(' Ground to excited state transition velocity quadrupole moments (Au):'):
              canprintlines = False
          if canprintlines:
              trs_mag_mom.append(line)

    trs_mag_mom=np.array(trs_mag_mom)    
    indi_mag_mom = []
    for i in range(1,len(trs_mag_mom)):
        indi_mag_mom.append(trs_mag_mom[i].split())

    indi_mag_mom = np.array(indi_mag_mom)
    
    test = []
    for i in range(len(indi_mag_mom)):
        test.append([indi_mag_mom[i][0], indi_mag_mom[i][1], indi_mag_mom[i][2],indi_mag_mom[i][3]])

    test = np.array(test)
    mag_mom_all = test[1:,1:4]
    return mag_mom_all


# Opening the gaussian file

try:
    ifile = sys.argv[1]          # opening the input file
    data = cclib.io.ccread(ifile)   
    
except IOError:
    print ("A valid gaussian output file must be provided")


# Obtaining the relevant quantities from the Gaussian output
    
osc_stregth = data.etoscs
natoms = data.natom
masses = data.atommasses
coord = data.atomcoords
dipole = data.moments

get_xyz(ifile)
#plot_dip_arrow(com_coord(coord,masses),gauss_trans_dip_mom(ifile)[42],ifile)
plot_dip_arrow(com_coord(coord,masses),gauss_trans_mag_mom(ifile)[0],ifile)

a=np.array(gauss_trans_mag_mom(ifile)[0],dtype=float)
b=np.array(gauss_trans_dip_mom(ifile)[0],dtype=float)

angle= np.arccos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)))
print (np.degrees(angle))

