#Protein/polmyer analysis functions

import MDAnalysis as mda
import MDAnalysis.analysis.distances as maa_dist
import statsmodels as stats
import math
import numpy as np
import pandas as pd
import sklearn


def get_protresd_list(prot_atoms, g2_atoms, dmax, universe):
    """Find all protein residues for which atoms that are within dmax. Here, contact is defined as any group two atoms
    that is within 4 Angstroms of the protein atom group (group 1). If that condition is satisfied, that residue is added to an array. 
    This calculation is only done on one frame.
    
    Inputs consisted of: 
    
    - prot_atoms: Atom group from MD Analysis containing the protein atoms 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups 
    - universe: MD Analysis variable containing information for the enitre trajectory
    
    Output consists of: 
    
    - mk: numpy array containing the residues that did contact group atoms"""
    
    ro = len(prot_atoms) # number of protein residues  
    
    cl = len(g2_atoms) # numbers of atoms in group two  
    
    dij_tri = np.zeros(shape=(ro,cl)) # Initialize a matrix the whose dimensions are (ro, cl)
    
    # Use MD Analysis to find distance between the protein atom groups and group two groups 
    dij_tri = maa_dist.distance_array(prot_atoms.positions
                                      , g2_atoms.positions, box=universe.trajectory.ts.dimensions)
    
    # Get the indices of the atoms that meet the 4 angstroms cutoff 
    exp_prot_atoms = np.any(dij_tri <= dmax, axis=1)
    
    # Get the corresponding residues that contact group two atoms 
    mk = np.array(prot_atoms[exp_prot_atoms].residues)
    
    return mk


def aa_frmcount(prot_atoms, g2_atoms, dmax, universe, start, end):
    """This function will output a dictionary of AA protein residue number and its corresponding frame count and occupancy. 
    
    Inputs consisted of: 
    
    - prot_atoms: Atom group from MD Analysis containing the protein atoms 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups 
    - universe: MD Analysis variable containing information for the enitre trajectory
    - start: Start of the trajectory (first frame number)
    - end: End of the trajectory (last frame number)
    
    Output is consisted of: 
    
    - aa_dict: A dictionary whose keys are residues that contacted group two atoms within the trajectory and
               values are an array containing the frame count and occupancy. Here, occupancy is defined as the 
               frame count of each residue divided by the total number of frames in trajectory block. Frame count 
               is the number of frames that a group two atom was less than dmax distance from any protein atom. 
    """
    
    aa_dict = {} # Initialize a dictionary to store the output. 
    
    laa = np.zeros(shape=len(prot_atoms.residues)) # Initialize a numpy array filled with zeros whose dim is no. of protein residues 
    
    br = np.array(prot_atoms.residues) # Store protein residues as given from MD Analysis 
    
    # This for loop goes through the trajectory block
    for ts in universe.trajectory[start:end]: 
        
        count = 0  # start counter at zero 
        
        # Use the get_protresd_list function to get residues that contact the group two atoms 
        bsres = get_protresd_list(prot_atoms, g2_atoms, dmax, universe) 
        
        # If the output array is empty, then no residues were within 4 Angstroms of group two atoms 
        if bsres.size == 0: 
            pass
        elif bsres.size != 0: # If array is not empty
    
            count += 1 # add to counter
            
            # Go through each residue 
            for i in bsres.flat:
        
                res_ind = np.where(br == i)[0]  # Get index for stored residue array that matches the residues from output in get_protresd_list
            
                laa[res_ind[0]] = laa[res_ind[0]] + count # update frame count to the specific residure 
                
    fin_res = np.where(laa != 0)[0] # Get indexes for residues with frame counts that are not zero
    
    # Calculate occupancy for each residue and store in a dictionary 
    for i in fin_res.flat:
        aa_dict[str(list(prot_atoms.residues[i:i+1])[0])] = [laa[i], laa[i]/(end - start)]
    
    return aa_dict 


def grptwocnt_aa(prot_atoms, g2_atoms, dmax, universe):
    """ This function calculates the number of PLGA monomers within 4 A of a BSA AA residue and 
    return a numpy array containing the no. of PLGA monomers for each BSA AA residue at one single frame. Make sure you use a pdb file 
    for system description to get unique chain identifiers. This calculation is only done on one frame. Universe is used, in case 
    the periodic boundary is present. 
    
    GROMACS command: gmx trjconv -s topol.tpr -f traj_comp.xtc  -o new_conf.pdb -dump 0
 
    Inputs consist of: 
    
    - prot_atoms: Atom group from MD Analysis containing the protein atoms 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups
    - universe: MD Analysis variable containing information for the enitre trajectory
    
    Outputs consist of: 
    
    - np.array(co_grpaa): For each amino acid group type, the total number of group two residues are given in a numpy array. 
    
    - la: numpy array containing the no. of group two (PLGA) monomers that meet the distance for each protein residue.
          Some array values will be empty because some residues will not contact the protein. 
    
    """
    
    gg_dict = {} # Initialize a dictionary to store one of the outputs
    
    # Use MD Analysis distance array function to get matrix of dim (# of group two atoms, # of protein atoms) 
    plga_mon = maa_dist.distance_array(g2_atoms.positions, prot_atoms.positions, box=universe.trajectory.ts.dimensions)
    
    # Return the indices of the elements in plga_mon that are less than dmax and are non-zero. Output is a tuple
    plga_mon_resd = np.nonzero(plga_mon < dmax)
    
    # get polymer residue numbers for each atom in contact with the protein atom group
    a = np.array(list(g2_atoms[plga_mon_resd[0]].resids))
    
    # get polymer segid numbers for each atom in contact with the protein atom group
    b = np.array(list(g2_atoms[plga_mon_resd[0]].segids))
    
    # concatenate the two above arrays 
    new_gt = np.array([a,b])
    
    # indices where there is a change in values within the resids polymer matrix 
    res_nsgt = np.where(new_gt[0][:-1] != new_gt[0][1:])[0]

    # indices where there is a change in values within the resids protein matrix 
    res_prot_gt = np.where(prot_atoms[plga_mon_resd[1]].resids[:-1] != prot_atoms[plga_mon_resd[1]].resids[1:])[0]

    # Get the polymer residue indexes not included in the prot residues index list 
    mi_sgt = res_nsgt[np.where(np.isin(res_nsgt,res_prot_gt) == False )]
    
    # Add it to the protein residue index list and sort in ascending order 
    act_gt = np.sort(np.append(res_prot_gt,mi_sgt))

    # Initialize zeros array to store the no. of PLGA residues that contact each protein residue 
    la = np.zeros(shape=len(prot_atoms.residues))
    
    for i in act_gt.flat: 
    
        m_i = prot_atoms[plga_mon_resd[1]].resids[i]
        #print(m_i)

        la[m_i-1] = la[m_i-1] + 1
        
    # Find all residues that have more than zero PLGA monomers 
    fin_res = np.where(la != 0)[0]
    
    # Save number of PLGA residues per AA 
    for i in fin_res.flat:
        gg_dict[str(list(prot_atoms.residues[i:i+1])[0])] = la[i]

    # Grouping of residues in Smith et al  
    aromatic_res = ['PHE', 'TRP', 'TYR', 'HIS']
    hydrophobic_res = ['ALA', 'ILE', 'LEU', 'VAL', 'GLY', 'PRO','PHE', 'TRP','MET']
    polar_res = ['ASN', 'CYS', 'GLN', 'SER', 'THR','TYR']
    neg_res = ['ASP', 'GLU']
    pos_res = ['ARG', 'HIS', 'LYS']

    frac_res = [neg_res, pos_res, polar_res, hydrophobic_res, aromatic_res]
    
    # For each amino acid type in frac_res, this code chunk saves the PLGA residue count in a list and sums them together to 
    # to get a total number of PLGA residues within the trajectory for each AA group in frac_res
    co_grpaa = []

    for row in frac_res:
        fr_list = []
        for j in range(len(row)):
            for key, value in gg_dict.items():
                if row[j] in key:
                    fr_list.append(value)
        co_grpaa.append(sum(fr_list))
        
    return np.array(co_grpaa), la


def gtwo_trjcnt(prot_atoms, g2_atoms, dmax, universe, start, end):
    """This function calcuates the average number of PLGA mononers per BSA AA group and no of PLGA residues for each AA for each traj. block. 
        The grptwocnt_aa is used to calculate the above values for each frame. 
        
    Inputs consist of: 
    
    - prot_atoms: Atom group from MD Analysis containing the protein atoms 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups
    - universe: MD Analysis variable containing information for the enitre trajectory
    - start: Beginning of trajectory block (first frame number) 
    - end: End of trajectory block (last frame number)
    
    Outputs consist of: 
    
    - aa_dict: dicitionary containing mean and pop. standard deviation of the no. of polymer residues for each different AA type 
    - l_final: the average no. of PLGA residues for each residue for each traj length, shape is (no. of prot_atoms residues)
    
    """
    
    # Different amino acid types 
    sf_lb = ["Negative", "Positive", "Polar", "Hydrophobic", "Aromatic"]
    
    # Initialize empty dictionary 
    aa_dict = {}
    
    # Initialize matrix whose dim is (prot_res, length of traj. block)
    laa = np.zeros(shape=(len(prot_atoms.residues),end-start))

    # Initialize matrix to store amino acid types for each frame
    trj_aa = np.zeros(shape=(len(sf_lb), end-start))
    
    # Initialize array 
    l_final = np.zeros(shape=(len(prot_atoms.residues)))
    
    # Make sure trajectory is at the beginning 
    universe.trajectory[start]
    
    count = 0
    
    for ts in universe.trajectory[start:end]: 
    
        # fr_pres: For each amino acid group type, the total number of group two residues is given 
        # hh_matx: No. of group two (PLGA) monomers that meet the distance for each protein residue
        fr_pres, hh_matx = grptwocnt_aa(prot_atoms, g2_atoms, dmax, universe)
        
        # No. of polymer residues that were in contact with the protein at each frame 
        laa[:,count] = hh_matx
 
        # For each amino acid type, save total number of PLGA residues at each frame for each AA group
        trj_aa[:, count] = fr_pres
        
        count += 1
    
    # Calculating the average no. of PLGA residues for each residue for each traj length
    for i in range(len(l_final)):
        l_final[i] = np.mean(laa[i,:])

    # get mean and pop. standard deviation of the no. of polymer residues for each different AA type 
    for i in range(len(sf_lb)):
        aa_dict[sf_lb[i]] = [np.mean(trj_aa[i,:]), np.std(trj_aa[i,:])]
    
    return aa_dict, l_final


def frac_cont(frm_count_dict):
    """ This function specifically takes as input the dictionary output from aa_frmcount function and 
    outputs a dictionary, nlkts, where each key is an amino acid type and whose corresponding value is as follows: 
    
    - nlkts: dictionary output     
    # co_grpaa[i] (Total no. of contacts for each Amino acid type ) 
    # tp_cnt[i] (total no. of amino acids for each type of amino acid type) 
    # norm_list[i] (total no. of contacts normalized by the protein surface fraction)
    # cont_l[i]] (Normalized fraction of contacts)
   
    """
    
    a_a = ["GLY","ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","SER","THR","CYS","TYR","ASN","GLN","ASP"
               ,"GLU","LYS","ARG","HIS"]
    
    # Grouping of residues in Smith et al  
    aromatic_res = ['PHE', 'TRP', 'TYR', 'HIS']
    hydrophobic_res = ['ALA', 'ILE', 'LEU', 'VAL', 'GLY', 'PRO','PHE', 'TRP','MET']
    polar_res = ['ASN', 'CYS', 'GLN', 'SER', 'THR','TYR']
    neg_res = ['ASP', 'GLU']
    pos_res = ['ARG', 'HIS', 'LYS']

    frac_res = [neg_res, pos_res, polar_res, hydrophobic_res, aromatic_res]
    sf_lbl = ["Negative", "Positive", "Polar", "Hydrophobic", "Aromatic"]
    
    # For each amino acid type in frac_res, this code chunk saves the frame count in a list and sums them together to 
    # to get a total frame count within the trajectory for each AA group in frac_res
    co_grpaa = []

    for row in frac_res:
        fr_list = []
        for j in range(len(row)):
            for key, value in frm_count_dict.items():
                if row[j] in key:
                    fr_list.append(value[0])
        co_grpaa.append(sum(fr_list))
        
    # This chunk of code gets an AA count from the above list, in order 
    # to get a total number of residues that contact BSA
    cpl_l = []

    for i in range(len(a_a)):
        count = 0
        for key, value in frm_count_dict.items():
            if a_a[i] in key:
                count += 1
        cpl_l.append(a_a[i]+" "+str(count))   

    # For each AA type in frac_res, this code chunk saves the count for each AA within 4 Angstroms of a PLGA trimer 
    # in a list based on the order in frac_res, then sums the counts to get a total number of AA for each AA type 
    tp_cnt = []   
    
    for row in frac_res:
        nw_l = []
        for i in range(len(row)):
            for j in range(len(cpl_l)):
                if row[i] in cpl_l[j]:
                    nw_l.append(int(cpl_l[j][4:6]))
        tp_cnt.append(sum(nw_l))           
    
    # Get the total count of AA that are within 4 A of PLGA oligomer
    bsum = len(frm_count_dict.keys())
    
    # The code chunk normalized the frame count of each AA group type by the protein surface fraction 
    # of each amino acid type contacted by a polymer surrogate.
    norm_list = []
    for i in range(len(co_grpaa)):
        lk = tp_cnt[i]/bsum
        if lk == 0: 
            norm_list.append(co_grpaa[i])
        elif lk != 0:
            norm_list.append(co_grpaa[i]/(tp_cnt[i]/bsum))
            
    # This conde chunk calculates the fractional contact based on the normalized frame count 
    cont_l = []
    nsum = sum(norm_list)
    for i in range(len(norm_list)):
        cont_l.append(norm_list[i]/nsum)
    
    #Save values in a dictionary 
    # Legend: co_grpaa[i] (Total no. of contacts) 
    # tp_cnt[i] (total no. of amino acids for each type of contact) 
    # norm_list[i] (total no. of contacts normalized by the protein surface fraction)
    # cont_l[i]] (Normalized fraction of contacts)
    nlkts = {}
    for i in range(len(sf_lbl)):
        nlkts[sf_lbl[i]] = [co_grpaa[i], tp_cnt[i], norm_list[i], cont_l[i]]
        
    return nlkts


# I want a list of total fraction of contacts where length is determined by no. of blocks and a dictionary 
# of contact groups as keys and list of fractional contacts as values(length of list will be no. of blocks)
def bavg_frac_cnt(no_of_blks, prot_atoms, g2_atoms, dmax, universe, no_surf, begin, final):
    """
    This function takes as inputs: 
    
    Inputs consist of: 
    
    - no_of_blks: number of blocks
    - prot_atoms: Atom group from MD Analysis containing the protein atoms 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups 
    - universe: MD Analysis variable containing information for the enitre trajectory
    - no_surf: total number of protein surface residues 
    - begin: Beginning of trajectory block (first frame number)
    - final: End of trajectory block (last frame number)
    
    The grptwo_trjcnt, aa_frmcount, and frac_count functions are used to calculate outputs. 
    
    Outputs consist of: 
    
    - ot_dab (dictionary): fraction of contacts within each specified block for the AA groups in
        ["Negative", "Positive", "Polar", "Hydrophobic", "Aromatic", "total_frac"]
        
    - plga_dict (dictionary): average no of PLGA monomers and the std deviation for each block for each AA group type
    
    - aares_npa: the avg no of PLGA residues, (no_of_blks, len of protein residues) shape
    
    """

    # No. of frames in each block 
    n_size = (final - begin)/no_of_blks
    
    # Initialize list for storing total fraction of contacts for each block 
    frcb = []
    
    # Initialize distiionary for fraction of contacts for each AA type 
    ot_dab = {}
    
    # Make sure trajectory is at the beginning 
    universe.trajectory[begin]
    
    # Different AA group types 
    sf_lbl = ["Negative", "Positive", "Polar", "Hydrophobic", "Aromatic", "total_frac"]

    # Initialize matrix for AA type for each block 
    blk_nparr = np.zeros(shape=((len(sf_lbl)-1),no_of_blks))
    
    # Initialize array that will store the main AA group types fraction contacts 
    plga_nparr = np.zeros(shape=((len(sf_lbl)-1),no_of_blks), dtype=object)
    
    # Initialize matrix to store avg no. of PLGA residues within dmax for each residue in each traj. block  
    aares_npa = np.zeros(shape=(no_of_blks), dtype=object)
      
    for i in range(no_of_blks):
        
        tpl = []
 
        start = universe.trajectory.frame
        print(start)
    
        end = int(start + n_size)
        print(end)
        
        # A dictionary of AA protein residue number and its corresponding frame count and occupancy for each block 
        hn_bcks = aa_frmcount(prot_atoms, g2_atoms, dmax, universe, start, end)
        
        # a dictionary where each key is an amino acid type and whose corresponding value is as follows:    
        # co_grpaa[i] (Total no. of contacts for each Amino acid type ) 
        # tp_cnt[i] (total no. of amino acids for each type of amino acid type) 
        # norm_list[i] (total no. of contacts normalized by the protein surface fraction)
        # cont_l[i]] (Normalized fraction of contacts)
        ff_dict = frac_cont(hn_bcks)
        
        # the average number of PLGA mononers per protein AA group and no of PLGA residues for each AA for each trajectory block 
        nk_dict, res_mat = gtwo_trjcnt(prot_atoms, g2_atoms, dmax, universe, start, end)
        
        # total number of protein residues that contact polymer residues 
        lk = len(hn_bcks.keys())
        
        # total fraction of contacts is calculated here    
        frcb.append(lk/no_surf)

        # Normalized fractional contacts for each trajectory block 
        for key, value in ff_dict.items():
            tpl.append(value[3])
            
        count = 0
        for key, value in nk_dict.items():
            plga_nparr[count,i] = np.array(value)
            count += 1
            
        blk_nparr[:,i] = tpl  
        
        # (no_of_blks, len of protein residues) shape, for each block and for each AA residue, the avg no of PLGA residues is recorded
        aares_npa[i] = res_mat
        
	# Make sure trajectory is at the beginning of the next block
        universe.trajectory[end]
    
    # Save fractional contacts for each AA group type, each element in the value array corresponds to a block 
    # calculated value
    for i in range(len(sf_lbl)-1):
        ot_dab[sf_lbl[i]] = blk_nparr[i,:]
        
    # Save average no of PLGA monomers and the std deviation for each block for each AA group type into dictionary
    plga_dict = {}
    for i in range(len(sf_lbl)-1):
        plga_dict[sf_lbl[i]] = plga_nparr[i,:]
    
    # total fraction of contacts within the specified blocks    
    ot_dab[sf_lbl[5]] = np.array(frcb)   
    
    # output legend 
    # ot_dab (dictionary): fraction of contacts within each specified block for the AA groups
    # plga_dict (dictionary): average no of PLGA monomers and the std deviation for each block for each AA group type
    # aares_npa: the avg no of PLGA residues, (5,583) shape
    
    return ot_dab, plga_dict, aares_npa


def prot_poly_cntmovie(prot_atoms, g2_atoms, dmax, universe, start, end):
    """This function calculates the contact matrix between protein AA residues and group 2 residues at each timestep 
    and returns a multi-dim numpy array saves contact map information for the entire trajectory block. A contact is saved as 
    a boolean value (0 or 1)
    
    Inputs consist of: 
    
    - prot_atoms: Atom group from MD Analysis containing the protein atoms 
    - g2_atoms: Atom group from MD Analysis for the second group 
    - dmax: Distance cutoff in Angstroms, used to define contact between the two groups
    - universe: MD Analysis variable containing information for the enitre trajectory
    - start: Beginning of trajectory block in picoseconds 
    - end: End of trajectory block in picoseconds 
    
    Output consists of: 
    
    - pp_mat: 3D numpy array (No. of protein residues, No. of group 2 residues, end-start) containing 
    contact map information for the entire trajectory block. A contact is saved as a boolean value (0 or 1)
    
    """
    
    # Make sure trajectory is at the beginning 
    universe.trajectory[start]
    
    cnt_un = 0 
    
    # Initialize output 3D array 
    pp_mat = np.zeros(shape=(end-start), dtype=object)
    
    
    for ts in universe.trajectory[start:end]:
    
        # Get distance matrix at each frame
        ro = len(prot_atoms)
        cl = len(g2_atoms)
        dij_tri = np.zeros(shape=(ro,cl))
        dij_tri = maa_dist.distance_array(g2_atoms.positions, prot_atoms.positions, box=universe.trajectory.ts.dimensions)

        # Initialize matrix to store contact info at each frame 
        matfr = np.zeros(shape=(len(prot_atoms.residues),len(g2_atoms.residues)))
    
        # Return the indices of the elements in pr_pol that are less than dmax and are non-zero. Output is a tuple
        pr_pol = np.nonzero(dij_tri < dmax)
    
        # Get residue IDs based on polymer atoms that contact protein 
        a = np.array(list(g2_atoms[pr_pol[0]].resids))
        
        # Get segment IDs based on polymer atoms that contact protein
        b = np.array(list(g2_atoms[pr_pol[0]].segids))
        
        # Concatenate the two arrays together 
        new_l = np.array([a,b])
    
        # indices where there is a change in values within the resids polymer matrix 
        res_nos = np.where(new_l[0][:-1] != new_l[0][1:])[0]

        # indices where there is a change in values within the resids protein matrix 
        res_prot = np.where(prot_atoms[pr_pol[1]].resids[:-1] != prot_atoms[pr_pol[1]].resids[1:])[0]

        # Get the polymer residue indexes not included in the prot residues index list 
        mi_s = res_nos[np.where(np.isin(res_nos,res_prot) == False )]
    
        # Add it to the protein residue index list and sort in ascending order 
        act = np.sort(np.append(res_prot,mi_s))

        for i in act.flat: 
    
            # Get protein residue index
            m_i = prot_atoms[pr_pol[1]].resids[i]
            
            # Get corresponding polymer residue index
            r_id1 = new_l[0][i]
    
            # Get corresponding polymer segment index
            r_seg1 = new_l[1][i]

            # Ensure the right chain contact is accounted for. 
            ar_res = np.array([20, 40])
    
            # Chain identfiers
            un_polys = np.unique(g2_atoms.segids)
    
            # Get polymer residue and segment index 
            if un_polys[0] == r_seg1:
                m_j = int(r_id1)
            elif un_polys[1] == r_seg1: 
                m_j = ar_res[0] + int(r_id1)
            elif un_polys[2] == r_seg1:
                m_j = ar_res[1] + int(r_id1)
        
            matfr[m_i-1,m_j-1] = 1
        
        # Save matrix for this frame 
        pp_mat[cnt_un] = matfr
    
        # update counter 
        cnt_un += 1
        
    return pp_mat


def AA_list_org(lorg_list):
    
    """List elements need have 'GLY  XX' as string format, where XX reps the number of GLY residues. Output is a
    sorted list of 'AA XX' according to the below order.  """
    
    hydrophobic_res = ['ALA', 'ILE', 'LEU', 'VAL', 'GLY', 'PRO','PHE', 'TRP','MET']
    polar_res = ['ASN', 'CYS', 'GLN', 'SER', 'THR','TYR']
    neg_res = ['ASP', 'GLU']
    pos_res = ['ARG', 'HIS', 'LYS']

    all_res = [pos_res, neg_res, polar_res, hydrophobic_res]
    #Change order of residues before making the bar graph
    # (1) Positively charged
    # (2) Negatively charged
    # (3) Polar residues 
    # (4) Hydrophobic residues 
    
    # This chunk of code sorts the counts of each AA that have 1001 or 1002 frame count based 
    # on the AA order in all_res
    arr_list = []

    for row in all_res:
        for i in range(len(lorg_list)):
            for j in range(len(row)):
                if row[j] == lorg_list[i][0:3]:
                    arr_list.append(lorg_list[i])
                    
    #This chunk of code splits the list arr_list to makes the AA: count of 1001 or 1002 frames data plottable 
    f_list = []
    fn_list = []
    for i in range(len(arr_list)):
        f_list.append(arr_list[i][0:3])
        fn_list.append(int(arr_list[i][5:]))
        
    return f_list, fn_list
