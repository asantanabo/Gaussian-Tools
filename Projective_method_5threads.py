import numpy as np
import sys
from cclib.io import ccopen
import scipy as sp
from scipy import linalg
import time
import queue, threading

my_queue = queue.Queue()


def storeInQueue(f):
  def wrapper(*args):
    my_queue.put(f(*args))
  return wrapper


@storeInQueue
def read_file(filn, no_mols):
    # parse file
    log = ccopen(filn)
    # features we want to extract:
    # Open the log files
    # Parse the relevant data
    mol = log.parse()
    # Get molecular orbitals. Need the transpose for our purposes.
    MOs=(mol.mocoeffs[0]).T

    # Get overlaps. These are symmetric so transpose not required in this case
    S=mol.aooverlaps
    # features we need depending on whether it's dimer or monomer
    if no_mols < 2:
        # Size of basis sets
        nbasis = mol.nbasis
        # Position of HOMO
        nhomo = mol.homos
        #print(nhomo)
        return nbasis, MOs, S, nhomo
    elif no_mols == 2:
        # Get eigenvalues of pair
        EvalsAB = mol.moenergies[0]
        Dpair=sp.linalg.cholesky(S)
        return MOs, Dpair, EvalsAB

@storeInQueue
def fill_MOs(Nbasis, nbasisA, MOsA, nbasisB, MOsB):
    MOs=np.zeros((Nbasis,Nbasis))
    # filling A and B into C local can be done in parallel not sequentially
    MOs[0:nbasisA,0:nbasisA]=MOsA
    MOs[nbasisA:Nbasis,nbasisA:Nbasis]=MOsB
    return MOs

@storeInQueue
def fill_S(Nbasis, nbasisA, SA, nbasisB, SB):
    S=np.zeros((Nbasis,Nbasis))
    # filling overlap matrix with A and B can be done in parallel not sequentially
    S[0:nbasisA,0:nbasisA]=SA
    S[nbasisA:Nbasis,nbasisA:Nbasis]=SB
    D = sp.linalg.cholesky(S)
    return D

# INPUT VARIABLES
# Read in molecule log files for projective method. Requires iop(3/33=1,6/7=3) in Gaussian header for calculation on each molecule + the pair
MOLA_proj=sys.argv[1]
MOLB_proj=sys.argv[2]
MOLAB_proj=sys.argv[3]
Degeneracy_HOMO=int(sys.argv[4])
Degeneracy_LUMO=int(sys.argv[5])    # =0 for non-degenerate, 2,3 etc for doubly, triply etc

if __name__ == "__main__":
    start_time = time.time()
    # start threads for reading in files in parallel
    t1 = threading.Thread(target=read_file, args=(MOLA_proj, 1))
    t2 = threading.Thread(target=read_file, args=(MOLB_proj, 1))
    t3 = threading.Thread(target=read_file, args=(MOLAB_proj, 2))

    t1.start()
    nbasisA, MOsA, SA, nhomoA = my_queue.get()
    #print("Mol A: All parameters have been extracted. ")
    t2.start()
    nbasisB, MOsB, SB, nhomoB = my_queue.get()
    #print("Mol B: All parameters have been extracted. ")
    t3.start()
    MOsAB, Dpair, EvalsAB = my_queue.get()
    t1.join()
    t2.join()
    Nbasis=nbasisA+nbasisB
    ############################################################
    # MOs = C local in  Kirkpatrick's method
    t4 = threading.Thread(target=fill_MOs, args=(Nbasis, nbasisA, MOsA, nbasisB, MOsB))
    t5 = threading.Thread(target=fill_S, args=(Nbasis, nbasisA, SA, nbasisB, SB))
    t4.start()
    MOs = my_queue.get()
    t5.start()
    D = my_queue.get()
    # Calculate upper diagonal matrix D, such that S=D.T*D for Loewdin orthogonalisation.
    # matrices NEEDED to get local basis sets
############################################################
    t4.join()
    t5.join()
    t3.join()


############################################################
# Orthogonalise MOs matrix and MOsAB matrix
# D unitary matrix, not normalised
    MOsorth=np.dot(D,MOs)

    MOspairorth=np.dot(Dpair,MOsAB)

    # Calculate the Fock matrix
    B=np.dot(MOsorth.T,MOspairorth)
    Evals=np.diagflat(EvalsAB)
    # inner np.dot: C pair Cloc).T
    F=np.dot(np.dot(B,Evals),B.T) # ranges from 0 to 1, already orthonormalised

############################################################
# Output the HOMO-HOMO and LUMO-LUMO coupling elements fromt the Fock matrix

if Degeneracy_HOMO==0:
    print("HOMO-HOMO coupling: ", F[nhomoB+nbasisA,nhomoA])

if Degeneracy_LUMO==0:
    print("LUMO-LUMO coupling: ", F[nhomoB+nbasisA+1,nhomoA+1])

# Degeneracies

if Degeneracy_HOMO!=0:
    F_deg_HOMO=F[nhomoB+nbasisA-Degeneracy_HOMO+1:nhomoB+nbasisA+1,nhomoA-Degeneracy_HOMO:nhomoA]
    print("HOMO-HOMO coupling", (np.sum(np.absolute(F_deg_HOMO**2)))/Degeneracy_HOMO**2)

if Degeneracy_LUMO!=0:
    F_deg_LUMO=F[nhomoB+nbasisA:nhomoB+nbasisA+Degeneracy_LUMO,nhomoA:nhomoA+Degeneracy_LUMO]
    print("LUMO-LUMO coupling", (np.sum(np.absolute(F_deg_LUMO**2)))/Degeneracy_LUMO**2)

print("Parallel run time: {} seconds".format(time.time() - start_time))
