'''
Code was adapted from
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb
'''


'''
---------------------
Import modules
---------------------
'''

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import os
import subprocess
import shutil
import signal
from contextlib import contextmanager

from Bio.PDB import PDBParser

from alphafold.relax import relax
from alphafold.relax import utils
from alphafold.common import protein, residue_constants

'''
---------------------
Define modified residues
---------------------
'''

MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

'''
---------------------
Define functions
---------------------
'''

def pdb_to_string(pdb_file, chains=None, models=[1]):
  '''read pdb file and return as string'''

  if chains is not None:
    if "," in chains: chains = chains.split(",")
    if not isinstance(chains,list): chains = [chains]
  if models is not None:
    if not isinstance(models,list): models = [models]

  modres = {**MODRES}
  lines = []
  seen = []
  model = 1
  for line in open(pdb_file,"rb"):
    line = line.decode("utf-8","ignore").rstrip()
    if line[:5] == "MODEL":
      model = int(line[5:])
    if models is None or model in models:
      if line[:6] == "MODRES":
        k = line[12:15]
        v = line[24:27]
        if k not in modres and v in residue_constants.restype_3to1:
          modres[k] = v
      if line[:6] == "HETATM":
        k = line[17:20]
        if k in modres:
          line = "ATOM  "+line[6:17]+modres[k]+line[20:]
      if line[:4] == "ATOM":
        chain = line[21:22]
        if chains is None or chain in chains:
          atom = line[12:12+4].strip()
          resi = line[17:17+3]
          resn = line[22:22+5].strip()
          if resn[-1].isalpha(): # alternative atom
            resn = resn[:-1]
            line = line[:26]+" "+line[27:]
          key = f"{model}_{chain}_{resn}_{resi}_{atom}"
          if key not in seen: # skip alternative placements
            lines.append(line)
            seen.append(key)
      if line[:5] == "MODEL" or line[:3] == "TER" or line[:6] == "ENDMDL":
        lines.append(line)
  return "\n".join(lines)

def relax_me(pdb_in, pdb_out, max_iterations, tolerance, stiffness, use_gpu):
  pdb_str = pdb_to_string(pdb_in)
  protein_obj = protein.from_pdb_string(pdb_str)
  amber_relaxer = relax.AmberRelaxation(
    max_iterations=max_iterations,
    tolerance=tolerance,
    stiffness=stiffness,
    exclude_residues=[],
    max_outer_iterations=3,
    use_gpu=use_gpu
  )
  relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=protein_obj)
  with open(pdb_out, 'w') as f:
      f.write(relaxed_pdb_lines)

# Get the length of the protein from the PDB file
def get_protein_length(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # Assuming we're interested in the first model
    protein_length = len(list(model.get_residues()))
    return protein_length

# Timeout code from https://stackoverflow.com/a/22348885
class TimeoutException(Exception): pass
@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

'''
---------------------
Relaxation parameters
---------------------
'''

max_iterations = 2000 #@param {type:"slider", min:0, max:5000, step:250}
#@markdown - Low values might never converge and crash and high values might also run for a long time. `0` = infinity is the AF2 default, however it might run for a long time.
tolerance = 2.39 # @param {type:"number"}
#@markdown - kcal/mol, the energy tolerance of L-BFGS.
stiffness = 10.0 # @param {type:"number"}
#@markdown - kcal/mol A^2, spring constant of heavy atom restraining potential.
#@markdown - **Note:** Descriptions and default values are taken AF2.
use_gpu = True # @param {type:"boolean"}
#@markdown - set "Runtime->Change runtime type" to a GPU in Google Colab if you activate this.

'''
---------------------
File management
---------------------
'''

# Get the path to a PDB file
try:
  unrelaxed_pdb_file = sys.argv[1]
except IndexError:
  print("Please provide the path to a PDB file.")
  sys.exit(1)

# Get the name of the PDB file and the folder it's in
pdb_file = unrelaxed_pdb_file.split("/")[-1]
unrelaxed_pdb_folder = unrelaxed_pdb_file.split("/")[-2]

# Get the time limit in seconds
try:
    time_limit_in_seconds = int(sys.argv[2])
except IndexError:
    time_limit_in_seconds = 60
    #print("Please provide the time limit in seconds.")
    #sys.exit(1)

# Check if the PDB file exists
if not os.path.isfile(unrelaxed_pdb_file):
  print(f"The PDB file {unrelaxed_pdb_file} does not exist.")
  sys.exit(1)
if not unrelaxed_pdb_file.endswith(".pdb"):
    print(f"The file {unrelaxed_pdb_file} is not a PDB file.")
    sys.exit(1)

# Make a folder for relaxed PDB files if it doesn't exist
relaxed_pdb_folder = unrelaxed_pdb_folder.split("/")[-1] + "_relaxed"
os.makedirs(relaxed_pdb_folder, exist_ok=True)

# Make a folder for unrelaxed PDB files that finished relaxation
unrelaxed_pdb_finished_folder = unrelaxed_pdb_folder + "_finished"
os.makedirs(unrelaxed_pdb_finished_folder, exist_ok=True)

# Make a folder for PDB files that failed to relax if it doesn't exist
failed_pdb_folder = unrelaxed_pdb_folder.split("/")[-1] + "_failed"
os.makedirs(failed_pdb_folder, exist_ok=True)

log_file = os.path.join(unrelaxed_pdb_folder, "relaxation_log.tsv")

'''
---------------------
Relax the protein
---------------------
'''

pdb_in = os.path.join(unrelaxed_pdb_folder, pdb_file)
relaxed_pdb_file = pdb_file[:-4] + "_relaxed.pdb"
pdb_out = os.path.join(relaxed_pdb_folder, relaxed_pdb_file)

protein_length = get_protein_length(pdb_in)

# Set a 60-second alarm for the relaxation process
try:
    with time_limit(time_limit_in_seconds):
            
            # Relax the protein
            relax_me(pdb_in=pdb_in,
                    pdb_out=pdb_out,
                    max_iterations=max_iterations,
                    tolerance=tolerance,
                    stiffness=stiffness,
                    use_gpu=use_gpu)
            
            # Print the status of the relaxation process and move the PDB file to the folder for finished PDB files
            shutil.move(pdb_in, os.path.join(unrelaxed_pdb_finished_folder, pdb_file))
            
            with open(log_file, "a") as f:
                f.write(f"{pdb_file}\tfinished\t{protein_length}\n")

# If the relaxation process takes too long, move the PDB file to the failed folder
except TimeoutException as e:
    
    # Print the status of the relaxation process and move the PDB file to the folder for failed PDB files
    shutil.move(pdb_in, os.path.join(failed_pdb_folder, pdb_file))
    
    with open(log_file, "a") as f:
        f.write(f"{pdb_file}\tfailed\t{protein_length}\n")
        