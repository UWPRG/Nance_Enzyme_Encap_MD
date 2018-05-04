import mdtraj as md
import os.path as op

# specify path to starting configuration file
data_dir = op.join(op.dirname(__file__), 'conf_data/')
conf_file = op.join(data_dir, 'PEG.pdb')

# unique residue types in conf_file
start_cap_name = 'mPEG'
repeat_name = 'PEG'
end_cap_name = 'PEGo'

# three reference atoms to translate and rotate monomer
# to polymerization
ref_atom_names = ['O', 'C1', 'C2']
# read conf_file with mdtraj
monomer = md.load(conf_file)
coordinates = monomer.xyz[0]
topology = monomer.topology

# select atom IDs for each unique segment
start_cap = topology.select(f'resname {start_cap_name}')
repeat = topology.select(f'resname {repeat_name}')
end_cap = topology.select(f'resname {end_cap_name}')

ref_atoms = [coordinates[topology.select('resname {} and name {}'
                                         ''.format(repeat_name, ref_atom))]
             for ref_atom in ref_atom_names]

# TODO: polymerize

