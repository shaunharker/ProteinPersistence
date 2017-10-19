# Set up environment
from Bio.PDB import *
import re
import ProteinPersistence
import json
import multiprocessing as mp
import sys

def Compute_PDB_PersistenceDiagrams(filename):
    # Load an example PDB file using BioPython
    p = PDBParser()
    structure = p.get_structure('',filename)
    # Create a list of entries [x, y, z, r] representing (x,y,z) coordinates of each atom
    # and van der Waals radius of that atom
    def atom_shortname(s):
        return s[re.search("[A-Z]", s).start()]

    def vanderWaalsRadius(atom_name):
        """ 
        Return van der Waals radius associated with 'atom_name'
        """
        return { "H" : 1.2, "C" : 1.7 , "N" : 1.55, "O" : 1.52, "S" : 1.8}[atom_name];
    xyzr_list = [ list(atom.get_coord()) + [vanderWaalsRadius(atom_shortname(atom.get_name()))] \
                  for atom in structure.get_atoms() ]
    # Compute 
    return ProteinPersistence.pdb2persistence(xyzr_list)

def ProcessFile(args):
    filename, path, savepath = args
    diagrams = Compute_PDB_PersistenceDiagrams( path + '/' + filename)
    with open(savepath + '/' + filename + '_persistence.json', 'w') as outfile:
        json.dump(diagrams, outfile)   

def PerformJobsWithMultiprocessing(manifest_filename, path, savepath):
    with open(manifest_filename) as f:
        manifest = f.readlines()
    pool = mp.Pool(processes=20)
    pool.map(ProcessFile, [ [filename, path, savepath] for manifest in manifest])

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: Compute_PDB_Diagrams.py <manifest_file> <path_of_pdb_files> <path_to_save_output>")
        exit(1)
    PerformJobsWithMultiprocessing(sys.argv[1], sys.argv[2], sys.argv[3])
