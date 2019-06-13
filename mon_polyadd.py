import os.path as op
from io import StringIO

import pandas as pd


data_dir = op.join(op.dirname(__file__), '/suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/PEG_chain/ffmaker_PEG')
ff_dir = op.join(op.dirname(__file__), '/suppscr/pfaendtner/cnyambura/NEE_home/BSA_Nano_Prep/polymer_force_field/conf_data/MOD_n_amber99sb-ildn.ff/')

# give file names for topology and force field
polymer_topolgy = op.join(data_dir, 'PEG_GMX.top')
force_field_nonbonded = op.join(ff_dir, 'ffnonbonded.itp')
force_field_bonded = op.join(ff_dir, 'ffbonded.itp')


cap2 = {'name': 'tPEG',
        'list': ['O3','C5','H10','H11','C6','H12','H13','O4','H14']
	}

monomer = {'name': 'PEG',
           'list': ['O2','C3','H6','H7','C4','H8','H9']
           }

cap0 = {'name': 'sPEG',
        'list': ['H1','O1','C1','H2','H3','C2','H4','H5']
        }

# read file into string for parsing
def parse_topology(filename, section, header):
    top_string = ''
    with open(filename) as top_file:
        for line in top_file:
            top_string += line
    # process section string for conversion to dataframe
    section_string = top_string.split(section)[1].split('\n\n')[0]
    pandas_readable = StringIO(section_string)
    # create section DataFrame
    section_df = pd.read_table(pandas_readable, sep='\s+',
                               comment=';', names=header)
    return section_df


# define headers for sections of interest
atoms_header = ['nr', 'type',  'resi',  'res',
                'atom', 'cgnr', 'charge', 'mass']

atom_types_header = ['name', 'bond_type', 'mass', 'charge',
                     'ptype', 'sigma', 'epsilon']

bonds_header = ['ai', 'aj', 'funct', 'r', 'k']

# get polymer topology DataFrames
atoms = parse_topology(polymer_topolgy, '[ atoms ]', atoms_header)
bonds = parse_topology(polymer_topolgy, '[ bonds ]', bonds_header)
atom_types = parse_topology(polymer_topolgy, '[ atomtypes ]',
                            atom_types_header)

for index, row in atoms.iterrows():

    if row['atom'] in monomer['list']:
        atoms.loc[index, 'restype'] = monomer['name']
    elif row['atom'] in cap0['list']:
        atoms.loc[index, 'restype'] = cap0['name']
    elif row['atom'] in cap2['list']:
        atoms.loc[index, 'restype'] = cap2['name']



# get force field topology DataFrames
ff_bonds = parse_topology(force_field_bonded, '[ bondtypes ]',
                          bonds_header)
ff_atom_types = parse_topology(force_field_nonbonded, '[ atomtypes ]',
                               atom_types_header)


for index, row in atom_types.iterrows():
    # attempt to assign atom type based on LJ parameters
    candidates = ff_atom_types.ix[
        (ff_atom_types['sigma'] == row['sigma'])
        & (ff_atom_types['epsilon'] == row['epsilon'])]

    # holler if there are no matches!
    if candidates.empty:
        print('Ooops! Could not find match for atom type {} in '
              'this forcefield! Add this atom manually.'.format(row['name']))
        continue
    # assign atom type if there is only match
    if len(candidates) == 1:
        atoms.loc[atoms['type'] == row['name'], 'ff_name'] = candidates['name'].values[0]
        continue


    atom_number = atoms[atoms['type'] == row['name']]['nr'].values[0]
    atom_names = str(
        atoms[atoms['type'] == row['name']]['atom'].values
    ).replace('[', '').replace(']', '')

    bond_check = bonds.ix[
        (bonds['ai'] == atom_number)
        | (bonds['aj'] == atom_number)]

    bond_candidates = ff_bonds.ix[
        (ff_bonds['r'] == bond_check['r'].values[0])
        & (ff_bonds['k'] == bond_check['k'].values[0])]

    if bond_candidates.empty:
        print("Found several candidate atomtypes based "
              "on LJ paramters, but no perfect bond mathches.\n")

        # quick last ditch assignment based on name
        if row['name'].upper() in candidates['name'].values:
            print("Approximate name match found. Would you like "
                  "to assign atomtype '{}' from the force field "
                  "to your atom(s) {}?".format(
                  row['name'].upper(), atom_names))
            choice = input("(y/n): ")
            if choice == 'y':
                atoms.loc[atoms['type'] == row['name'], 'ff_name'] = row['name'].upper()
                continue

        print("Here is a list of other candidates.")
        print("Please choose an atom name for '{}' "
              "from the following list:".format(row['name']))
        for idx, name in enumerate(candidates['name'].values):
            print('{}. {}'.format(idx, name))
        selection = int(input('number: '))
        atom_types.loc[index, 'ff_name'] = candidates['name'].values[selection]

# Alright.... let's say we've got atom types correctly assigned by this point.
# We'll have to do some manual editing, but we're close enough for government
# work. Now to format these things in the correct group.

for new_residue in atoms['restype'].unique():
    res_df = atoms.loc[atoms['restype'] == new_residue, :]
    header = '[ {} ]\n [ atoms ]\n'.format(new_residue)
    res_string = [header]
    for index, row in res_df.iterrows():
        start_index = res_df.index[0]
        name = row['atom']
        ff_type = row['ff_name']
        charge = row['charge']
        id = row['nr'] - start_index

        row_string = f"{name:>6}    {ff_type:<3}{charge:17.5f} {id:5}\n"
        res_string.append(row_string)

    res_string.append(' [ bonds ]\n')

    for index, row in bonds.iterrows():
        atom1 = atoms.loc[atoms['nr'] == row['ai'], 'atom'].values[0]
        atom2 = atoms.loc[atoms['nr'] == row['aj'], 'atom'].values[0]

        if (atom1 in res_df['atom'].values) \
                & (atom2 in res_df['atom'].values):

            row_string = f'{atom1:>6} {atom2:>5}\n'
            res_string.append(row_string)

        elif atom1 in res_df['atom'].values:
            atom2 = '+' + atom2
            row_string = f'{atom1:>6} {atom2:>5}\n'
            res_string.append(row_string)
        elif atom2 in res_df['atom'].values:
            atom1 = '-' + atom1
            row_string = f'{atom1:>6} {atom2:>5}\n'
            res_string.append(row_string)

    res_top = ''.join(res_string)
    outfile = "{}_top.itp".format(new_residue)
    full_out_path = op.join(data_dir, outfile)

    with open(full_out_path, "w") as text_file:
        print(res_top, file=text_file)

