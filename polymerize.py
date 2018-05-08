import math
import os.path as op

import numpy as np
import pandas as pd
import mdtraj as md

# specify path to starting configuration file
data_dir = op.join(op.dirname(__file__), 'conf_data/')
conf_file = op.join(data_dir, 'PEG.pdb')
outfile = 'polymer.pdb'

# unique residue types in conf_file
start_cap_name = 'mPEG'
repeat_name = 'PEG'
n = 5  # number of repeat monomers
flip_repeat = np.pi  # angle to rotate each repeat unit
end_cap_name = 'PEGo'

# three reference atoms to translate and rotate monomer
# to polymerization
ref_atom_names = ['resname mPEG and name O',
                  'resname PEG and name O',
                  'resname PEGo and name O']
# read conf_file with mdtraj
monomer = md.load(conf_file)
coordinates = monomer.xyz[0]
topology = monomer.topology
table, bonds = topology.to_dataframe()

# select atom IDs for each unique segment
start_cap = topology.select(f'resname {start_cap_name}')
repeat = topology.select(f'resname {repeat_name}')
end_cap = topology.select(f'resname {end_cap_name}')

ref_atoms = [topology.select(f'{ref_atom}')
             for ref_atom in ref_atom_names]

# POLYMERIZE
# 1. translate so the first ref atom is at the origin
coordinates -= coordinates[ref_atoms[0]]

# 2. pivot around ref_atom[0] so third ref atom is on the x axis z=0
def rotate(coords, axis, theta):
    """
    Return the new coordinates after counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = axis / np.linalg.norm(axis)
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rotation = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    return np.dot(coords, rotation.T)


# calculate normal vector and angle between vector defined by
# ref_atoms[0] and ref_atoms[2] and the target orientation
original = coordinates[ref_atoms[2]]
target = np.array([1, 0, 0])
normal = np.cross(original, target)[0]
angle = math.acos(np.dot(original, target)
                  / (np.linalg.norm(original) * np.linalg.norm(target)))
# rotate coordinates
coordinates = rotate(coordinates, normal, angle)

# 3. rotate so second ref atom has z=0
original = np.array([0,
                     coordinates[ref_atoms[1]][0, 1],
                     coordinates[ref_atoms[1]][0, 2]])
original = original.reshape(1, 3)
target = np.array([0, 1, 0])
normal = np.cross(original, target)[0]
angle = math.acos(np.dot(original, target)
                  / (np.linalg.norm(original) * np.linalg.norm(target)))
# rotate coordinates
coordinates = rotate(coordinates, normal, angle)

# add positional data to table
table['x'] = coordinates[:, 0] * 10
table['y'] = coordinates[:, 1] * 10
table['z'] = coordinates[:, 2] * 10

# list to hold residue dataframes
my_polymer = [table.iloc[start_cap]]

# get x and y values for repeat translation
x_shift, y_shift, _ = (coordinates[ref_atoms[2]]
                       - coordinates[ref_atoms[1]])[0] * 10
# 4. append repeats to my_polymer
for repeat_idx in range(0, n):
    indices = [i + repeat_idx * len(repeat)
               for i in repeat]
    # copy monomer info from original table
    new_monomer = pd.DataFrame(table.iloc[repeat].values,
                               index=indices,
                               columns=table.columns)

    # translate the coordinates, and flip if necessary
    new_monomer.x += x_shift * repeat_idx
    if flip_repeat and repeat_idx % 2 == 1:
        # translate in the y direction
        new_monomer.y += y_shift
        # flip around ref2 - ref0 axis
        axis = (coordinates[ref_atoms[2]]
                - coordinates[ref_atoms[0]])[0]
        flipped = rotate(new_monomer.loc[:, ['x', 'y', 'z']].values,
                         axis, flip_repeat)
        new_monomer[['x', 'y', 'z']] = flipped

    new_monomer.serial += repeat_idx * len(repeat)
    new_monomer.resSeq += repeat_idx
    my_polymer.append(new_monomer)

# 5. add capping group
indices = [i + (n - 1) * len(repeat)
           for i in end_cap]
# copy end_cap info from original table
new_cap = pd.DataFrame(table.iloc[end_cap].values,
                       index=indices,
                       columns=table.columns)

# translate the coordinates, and flip if necessary
new_cap.x += x_shift * (n - 1)
if n % 2 == 0 and flip_repeat:
    # translate by y distance from reference axis
    new_cap.y += y_shift
    # flip around ref2 - ref0 axis
    axis = (coordinates[ref_atoms[2]]
            - coordinates[ref_atoms[0]])[0]
    flipped = rotate(new_cap.loc[:, ['x', 'y', 'z']].values,
                     axis, flip_repeat)
    new_cap[['x', 'y', 'z']] = flipped

new_cap.serial += (n - 1) * len(repeat)
new_cap.resSeq += (n - 1)

# convert list to polymer DataFrame
my_polymer.append(new_cap)
polymer = pd.concat(my_polymer)

pdb_lines = [f'ATOM  {atom.serial:>5} {atom.name:<4} {atom.resName:>4} '
             f'{atom.resSeq:>4}  {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}'
             for atom in polymer.itertuples()]

pdb_string = '\n'.join(pdb_lines)
full_out_path = op.join(data_dir, outfile)
with open(full_out_path, "w") as text_file:
    print(pdb_string, file=text_file)
