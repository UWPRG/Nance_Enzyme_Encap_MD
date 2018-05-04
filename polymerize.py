import mdtraj as md
import os.path as op

# specify path to starting configuration file
data_dir = op.join(op.dirname(__file__), 'conf_data/')
conf_file = op.join(data_dir, 'PEG.pdb')

# read into mdtraj
monomer = md.load(conf_file)
