# write function to average across frames to give ensembled averaged cosine theta values at each N - 1 value
import MDAnalysis as mda
import statsmodels as stats
import math
import numpy as np 
import pandas
import sklearn
import scipy as scipy
import warnings 
# Use sklearn to do fitting
from sklearn.linear_model import LinearRegression
from scipy import stats
from scipy.optimize import leastsq

def pers_length(polymer_atoms, n_monomers):
    """ This function takes the polymer atoms and number of monomers and outputs the polymer-averaged cosine theta (dot product) values 
    at specific arc lengths (ds) in a numpy array. This function uses the center of mass points along the polymer chain to calc.
    the dot product of the vector from the 1st to 2nd residue and the vector from the 1st to 3rd, 1st to 4th and so on. 
    
    The structure of the output vec_poly is outlined below: 
    
    # Dimensions are (3 rows, N - 1 columns), N is n_monomers
    # row 1: Polymer-averaged cosine theta (dot product) values at each arc length value (length of vector is N - 1)
    # row 2: Evenly spaced integer values beginning from 1 to N - 1
    # row 3: Polymer-averaged angle values in degrees at each arc length value (Take the inverse cosine of each element in row 1)
    
    The polymer_atoms variable is an output from the atom selection class in MD Analysis. n_monomers
    variable is an integer value. """
    
    # Initialize a zeros matrix of 3 rows, N - 1 columns to store key values
    vec_poly = np.zeros(shape=(3,n_monomers-1), dtype=float)
    
    # Inititalize a evenly spaced values with length of N - 1
    len_vec = np.arange(n_monomers-1)

    # Store them in the vec_poly matrix
    vec_poly[1,:] = len_vec

    # Initialize a counter 
    count = 0 
    
    # Initialize a vector to store each cosine theta value, dtype is object because 
    # as the code moves down the polymer chain to calc the dot product at each arc length,
    # the reference starting point shift down the polymer chain to get more dot product samples
    sv_ply = np.zeros(shape=(n_monomers-1), dtype=object)
    
    # For loop that is the length of the monomers 
    for i in range(n_monomers):
        
        # Add 1 to counter
        count += 1
        
        # This array dim(N - count) will contain the dot product values for each round
        ds_cor = np.zeros(shape=(n_monomers-count))
        
        # For loop that is N - count monomers
        for j in range(n_monomers - count):
        
            # Inititalize a evenly spaced values with length of N - count
            jh = np.arange(n_monomers - count)
            
            # Add 1 to all values for the purpose of selecting the next residue for the dot product calculation
            jh += count+1

            # Use MD Analysis to select the first reference residue
            n6_mon1 = polymer_atoms.select_atoms("resid "+str(count))

            # Use MD Analysis to select the second reference residue
            n6_mon2 = polymer_atoms.select_atoms("resid "+str(jh[j]))

            # This if statement ensure calc. of the first dot product at the first arc length
            # since this calculation is repeated for the next first arc length/dot product calc
            # started after the second round of the initial first for loop 
            if j == 0: 
                
                # Calc. distance between center of masses 
                v1 = n6_mon1.center_of_mass() - n6_mon2.center_of_mass()
                    
                # normalize the distance to form normalized vector
                v1_norm = v1/(np.linalg.norm(v1))
            
                # take the dot product of the vetor with its self to get 
                # cosine theta value at the first arc length
                ds_cor[j] = v1_norm.dot(v1_norm)

             # This if statement calc the dot product at all N - count arc length values
            elif j != 0:
                
                # Calc. distance between center of masses
                v2 = n6_mon1.center_of_mass() - n6_mon2.center_of_mass()
                
                # Calc. distance between center of masses
                v2_norm = v2/(np.linalg.norm(v2))
                
                # take the dot product of the first arc length vector with the next to get 
                # cosine theta value at each subsquent arc length
                ds_cor[j] = np.dot(v1_norm, v2_norm)

            # Store the saved dot product values in the sv_ply vector
            sv_ply[i] = ds_cor
    
    # This vector will contain the polymer averaged dot product values 
    cor_avg = []
 
    # There are N - 1 cosine theta values
    # this loop will store each dot product sample in a list, take the mean and store it in cor_avg
    for j in range(n_monomers-1):
        lss = []
        for i in sv_ply.flat:
            try:
                lss.append(i[j])
            except IndexError:
                pass
        # apeend average cosine theta value to list     
        cor_avg.append(np.mean(lss))
    
    # Turn cor_avg into a numpy array 
    nm = np.array(cor_avg)

    # This vector will contain the polymer averaged angles at each ar length value
    # This is simply calculated by taking the inverse cosine of the dot product values 
    ang_vg = []
    
    # This loop does the inverse cosine calc, floating point cutoff works well but could be better
    for i in nm.flat:
        if i >= float(0.99):
            ang_vg.append(0)
        elif i <= float(0.99):
            ang_vg.append(math.degrees(math.acos(i)))

    # Store Polymer-averaged cosine theta (dot product) values
    vec_poly[0,:] = nm 
    
    # Store Polymer-averaged angle values in degrees
    vec_poly[2,:] = np.array(ang_vg)
    
    return vec_poly




def mean_sq_e2e(polymer_atoms, universe, n_monomers, start, end):
    
    # Initialize a row vector to store the end to end distance vector at each frame 
    e2e_ens = np.zeros(shape=(end-start), dtype=object)

    # Magnitudes (or distances) end to end distance vector to be stored here
    # row 1: squared end to end distance at each frame 
    # row 2: End to end distance at each frame 
    dis_e2e = np.zeros(shape=(2, end-start))
    
    count_e2e = 0 
    
    c_res = 0
    
    for ts in universe.trajectory[start:end]:
        
        # Get the position vector for COM of the first residue in the polymer atom group
        pres_start = polymer_atoms.select_atoms("resid "+str(c_res+1)).center_of_mass()
        
        # Get the position vector for COM of the last residue in the polymer atom group
        pres_end = polymer_atoms.select_atoms("resid "+str(n_monomers)).center_of_mass()
        
        # Calc end to end distance vector
        e2e = pres_end - pres_start
        
        # Square each distance vector so (x^2 + y^2 + z^2) to get squared e2e distance at current frame
        dis_e2e[0,count_e2e] = np.linalg.norm(e2e)**2
        
        # Take square root of squared e2e distance to get e2e distance at current frame
        dis_e2e[1,count_e2e] = np.linalg.norm(e2e)
        
        # Store e2e vector at current frame
        e2e_ens[count_e2e] = e2e
        
        count_e2e += 1

    # Outputs: 
    # e2e_ens: squared end to end distance vector at each frame 
    # dis_e2e: end to end distance (in Angstroms) and mean squared e2e at each frame 
    return e2e_ens, dis_e2e



def hydro_rad_poly(polymer_atoms, universe, n_monomers, start, end):
    
    # Initialize a zeros matrix of 2 rows, N - 1 columns to store key values
    # Row 1: Rh at each frame 
    rh_frame = np.zeros(shape=(end-start))
    
    #store sums
    sum_rhstore = np.zeros(shape=(n_monomers-1))
    
    cnn = 0
    
    for ts in universe.trajectory[start:end]:
        
        count = 0
    
        # For loop that is the length of the monomers 
        for i in range(n_monomers-1):
        
            # Add 1 to counter
            count += 1
            
            #print(count)
        
            # This array dim(N - count) will contain the dot product values for each round
            r_vec = np.zeros(shape=(n_monomers-count))
        
            # For loop that is N - count monomers
            for j in range(n_monomers - count):
        
                # Inititalize a evenly spaced values with length of N - count
                jh = np.arange(n_monomers - count)
            
                # Add 1 to all values for the purpose of selecting the next residue for the dot product calculation
                jh += count+1

                # Use MD Analysis to select the first reference residue
                n6_mon1 = polymer_atoms.select_atoms("resid "+str(count))

                #print(count)
                
                # Use MD Analysis to select the second reference residue
                n6_mon2 = polymer_atoms.select_atoms("resid "+str(jh[j]))
                
                #print(jh)

                # Calc. distance between center of masses
                v1 = n6_mon1.center_of_mass() - n6_mon2.center_of_mass()
                
                # store the inverse of these distances, same as mda.distance.distance_array
                r_vec[j] = 1/(np.linalg.norm(v1))
                
                #print(r_vec)
            
            sum_rhstore[count-1] = np.sum(r_vec)
        
        #print(sum_rhstore)    
            
        #rh_frame[0,cnn] = (1/((n_monomers)**2))*(np.sum(sum_rhstore))
        
        rh_frame[cnn] = 1/((1/((n_monomers)**2))*(np.sum(sum_rhstore)))

        cnn += 1
        
    return rh_frame




def rh_block_avg(no_of_blks, polymer_atoms, universe, begin, final):
    
     # Block size 
    n_size = int((final - begin)/no_of_blks)
    
    # Initialize dictionary 
    ot_dab = {}
    
    # Shift the start of the trajectory to the begin variable 
    universe.trajectory[begin]
    
    # Using MD Analysis, I can get the number of polymer monomers in the polymer_atoms input variable. 
    n_monomers = len(np.unique(polymer_atoms.resids))
    
    inv_sto = np.zeros(shape=(no_of_blks), dtype=object)

    for nb in range(no_of_blks):
        
        # Shift the start of the trajectory to the start variable
        start = universe.trajectory.frame
        print(start)
    
        # Define end of block 
        end = int(start + n_size)
        print(end)
        
        # Initialize array to store Rh for each block
        inv_frame = np.zeros(shape=(end-start), dtype=object)
    
        cnn = 0
    
        for ts in universe.trajectory[start:end]:
        
            sum_rhstore = np.zeros(shape=(n_monomers-1), dtype=object)
            #print(sum_rhstore.shape)
            
            co_t = 0
    
            # For loop that is the length of the monomers 
            for i in range(n_monomers-1):
        
                # Add 1 to counter
                co_t += 1
        
                # This array dim(N - count) will contain the dot product values for each round
                r_vec = np.zeros(shape=(n_monomers-co_t))
        
                # For loop that is N - count monomers
                for j in range(n_monomers - co_t):
        
                    # Inititalize a evenly spaced values with length of N - count
                    jh = np.arange(n_monomers - co_t)
            
                    # Add 1 to all values for the purpose of selecting the next residue for the dot product calculation
                    jh += co_t+1

                    # Use MD Analysis to select the first reference residue
                    n6_mon1 = polymer_atoms.select_atoms("resid "+str(co_t))

                    #print(polymer_atoms.select_atoms("resid "+str(co_t)).resids)
                
                    # Use MD Analysis to select the second reference residue
                    n6_mon2 = polymer_atoms.select_atoms("resid "+str(jh[j]))
                    
                    #print(polymer_atoms.select_atoms("resid "+str(jh[j])).resids)

                    # Calc. distance between center of masses 
                    v1 = n6_mon1.center_of_mass() - n6_mon2.center_of_mass()
                
                    # store the inverse of these distances
                    # This gives the same values as mda.distance.distance_array
                    r_vec[j] = 1/(np.linalg.norm(v1))
                
                    #print("r_vec ="+str(r_vec))
            
                sum_rhstore[co_t-1] = r_vec
                
                #print("sum_rhstore ="+str(sum_rhstore))
                
            inv_frame[cnn] = sum_rhstore
            
            cnn += 1
            
        # Enseble averaging of inverse distances prior to summation
        conn = 0
        s_rh = np.zeros(shape=(n_monomers-1), dtype=object)

        for ns in range(n_monomers-1):
            conn += 1
            ln = np.zeros(shape=(n_size,n_monomers-conn))

            rvc = np.zeros(shape=(n_monomers-conn))
            
            #print("ns = "+str(ns))
    
            for kl in range(n_size):
                #print("kl = "+str(kl))
                ln[kl] = inv_frame[kl][ns]
        
            #print(ln)
            for k in range(n_monomers-conn):
                #print(ln[:,k])
                #print(np.mean(ln[:,k]))
                rvc[k] = np.mean(ln[:,k])
        
            #print(rvc)
  
            s_rh[ns] = np.sum(rvc)
    
        #print(s_rh)

        # Calc time averaged Rh
        inv_sto[nb] = 1/((1/(n_monomers**2))*np.sum(s_rh))
        
        # Shift the start to the next trajectory block 
        universe.trajectory[end]

    return inv_sto



def orientation_order_param(polymer_atoms, universe, n_monomers, start, end):
    
    # Initialize a zeros matrix of 2 rows, N - 1 columns to store key values
    # Row 1: Rh at each frame 
    oo_frame = np.zeros(shape=(end-start))
    
    c_res = 0
    
    c_n2 = 0 

    for ts in universe.trajectory[start:end]:
    
        # Get the position vector for COM of the first residue in the polymer atom group
        pres_start = polymer_atoms.select_atoms("resid "+str(c_res+1)).center_of_mass()
        
        #print(polymer_atoms.select_atoms("resid "+str(c_res+1)).resids)
        
        # Get the position vector for COM of the last residue in the polymer atom group
        pres_end = polymer_atoms.select_atoms("resid "+str(n_monomers)).center_of_mass()
        
        #print(polymer_atoms.select_atoms("resid "+str(n_monomers)).resids)
        
        # Calc end to end distance vector
        e2e = pres_end - pres_start
        
        e2e_norm = e2e/(np.linalg.norm(e2e))
        
        # Inititalize a evenly spaced values with length of N - count
        jh = np.arange(n_monomers-1)
        
        jh += 2
        
        #print(jh)
        
        cosine_vals = np.zeros(shape=(n_monomers-1))
    
        # For loop that is the length of the monomers 
        for i in range(len(jh)):
            
            poly_mon = polymer_atoms.select_atoms("resid "+str(jh[i])).center_of_mass()
            
            if i == 0:
                
                oo_vec = poly_mon - pres_start
                
                oov_norm = oo_vec/(np.linalg.norm(oo_vec))
                
                # I get negative values, if I don't take absolute values 
                cosine_vals[i] = np.absolute(((np.dot(e2e_norm, oov_norm))**2) - 0.5)
                
            elif i != 0:
                                  
                poly_fmon = polymer_atoms.select_atoms("resid "+str(jh[i]-1)).center_of_mass()
                                  
                poly_smon = polymer_atoms.select_atoms("resid "+str(jh[i])).center_of_mass()
                                 
                oo_vec = poly_smon - poly_fmon
                
                oov_norm = oo_vec/(np.linalg.norm(oo_vec))
                
                # I get negative values, if I don't take absolute values
                cosine_vals[i] = np.absolute(((np.dot(e2e_norm, oov_norm))**2) - 0.5)
                
                #print(cosine_vals)
                
        oo_frame[c_n2] = (3/(2*(n_monomers-1)))*np.sum(cosine_vals)
        
        #print(oo_frame)
        
        c_n2 += 1
    
    return oo_frame


def get_rg_pers_poly(polymer_atoms, universe, start, end):
    """This function takes as inputs: 
    
    - polymer_atoms variable, an output from the atom selection class in MD Analysis.
    - the MD Analysis universe variable that contains trajectory information (Make sure you have removed PBC and made molecule whole)
    - Start time of the trajectory in picoseconds 
    - End time of the trajctory in picoseconds
    
    Once received, this function uses the pers_length function to retrieve polymer averaged cosine theta values at their 
    specified arc length values for each frame and stored in a output matrix, corr_v, since polymer atom coordinates shift at each frame.
    
    Those cosine theta values at each frame are averaged for each arc length value to get the time averaged, polymer averaged cosine theta 
    values. Those values are also used to calculate the standard deviation of the dot products for each arc length value. Both vectors are 
    stored in a output matrix, v_poly
    
    Radius of gyration calculation at each frame is also performed using the MD Analysis Rg function and stored in the
    output matrix, rg_ens. The time averaged radius of gyration is a floating point value stored in the output variable, avg_rg
    
    The structure of the outputs are given below:
    
    rg_ens:
            - Dim is 1 row, end - start columns, where each element is the radius of gyration of polymer_atoms at each frame
            
    v_poly:
            - Dim is 4 rows, N - 1 columns
              row 1: Time averaged, polymer averaged cosine theta values at each arc length value 
              row 2: Pop. standard deviation of each time averaged cosine theta at each arc length value 
              row 3: Time averaged, polymer averaged angle values in degrees at each arc length value 
              row 4: Evenly spaced integer values beginning from 1 to N - 1
    
    corr_v: 
            - Dim is N - 1 rows, end - start columns, where each row corresponds to each arc length value
              and each polymer averaged cosine theta values at each frame are stored in the columns. 
    
    avg_rg:
            - Integer value, time averaged radius of gyration
    """
    
    # Using MD Analysis, I can get the number of polymer monomers in the polymer_atoms input variable. 
    n_monomers = len(np.unique(polymer_atoms.resids))
    
    # Initialize a row vector to store the radius of gyration and squared radius of gyration at each frame 
    rg_sq_ens = np.zeros(shape=(2, end-start))
    
    # Initialize matrix to store polymer averaged cosine theta values at each frame
    corr_v = np.zeros(shape=(n_monomers-1,end-start))
    
    # Initialize matrix to store polymer averaged angle values (in degrees) at each frame
    angle_v = np.zeros(shape=(n_monomers-1,end-start))
    
    # Initialize matrix to store time averaged cosine theta for each arc length value
    v_poly = np.zeros(shape=(4,n_monomers-1))
    
    # initialize counter
    count_rg = 0
    
    # Move trajectory start time to start variable 
    universe.trajectory[start]
    
    for ts in universe.trajectory[start:end]:
        
        # Get polymer-averaged cosine theta values at this frame 
        # Even though the pers_length function cannot be fed the universe variable, the polymer
        # atom coordinates will update for the polymer atoms because of MD Analysis selection class
        p_mat = pers_length(polymer_atoms, n_monomers)

        # store polymer-averaged cosine theta values at this frame
        corr_v[:,count_rg] = p_mat[0]
        
        # store polymer averaged angle values at this frame
        angle_v[:,count_rg] = p_mat[2]
        
        # the radius of gyration at this frame 
        rg_sq_ens[0, count_rg] = polymer_atoms.radius_of_gyration()
        
        # the squared radius of gyration at this frame
        rg_sq_ens[1, count_rg] = (polymer_atoms.radius_of_gyration())**2
        
        # update counter 
        count_rg += 1 
        
    # store evenly spaced integer values beginning from 1 to N - 1 from pers_length function
    v_poly[3,:] = p_mat[1]
    
    # For loop calculates time averaged values 
    for i in range(n_monomers-1):
        
        # Calc time averaged cosine value at each arc length
        v_poly[0,i] = np.mean(corr_v[i,:])
        
        # Calc time std dev for each cosine theta are each arc length 
        v_poly[1,i] = np.std(corr_v[i,:])
        
        # Calc time averaged angle value (in degrees) at each arc length
        v_poly[2,i] = np.mean(angle_v[i,:])
    
    # Time averaged radius of gyration 
    avg_rg = np.sqrt(np.mean(rg_sq_ens[1]))
    
    return  rg_sq_ens, v_poly, corr_v, avg_rg 

def bavg_pers_cnt(no_of_blks, polymer_atoms, universe, len_bnd, fit_pnts, begin, final):
    """This function takes as inputs:

    - no_of_blks variable that defines the number of blocks
    - polymer_atoms variable, an output from the atom selection class in MD Analysis.
    - the MD Analysis universe variable that contains trajectory information (Make sure you have removed PBC and made molecule whole)
    - len_bnd variable that defines the length (in Angstroms) of the arc length 
    - fit_pnts dictates the number of points to fit, in order to calculate the persistence length. The number of points used to fit a 
      line to ln(cos(theta)) vs. arc length values dictated by the user will be the same for all blocks. 
    - Start time of the trajectory in picoseconds 
    - End time of the trajctory in picoseconds
    
    Once received, this function uses the get_rg_pers_poly function to retrieve the average radius of gyration and 
    the time averaged cosine theta values at each arc length value for each trajectory block. The time averaged values are then 
    used to fit to the function that relates the persistence length and the flexural rigidity(cosine theta). A linear fit on 
    ln(cos(theta)) vs. arc length, along with user defined number of points for fitting(fit_pnts), will provide the persistence length (Lp)
    values for each block. 
    
    As the functions is running, the start and end of the block, Lp, error in Lp from the fit and pearson coefficient are printed.  
    
    The structure of the outputs are given below:
    
    ot_dab:
            - Dictionary (2 keys, no_of_blks values per key), contains the radius of gyration and Lp at for each block. 
            
    mod_res:
            - Dim is 4 rows, no_of_blks columns
              row 1: Lp for each trajectory block  
              row 2: Error in Lp from fit [Angstroms], 95% Confidence Interval for each trajectory block 
              row 3: Model slope value from fit for each trajectory block 
              row 4: Mean squared error in Lp fit for each trajectory block     
    
    """

    # Website reference: https://www.statisticshowto.datasciencecentral.com/
    # probability-and-statistics/descriptive-statistics/sample-variance/
    
    # Block size 
    n_size = (final - begin)/no_of_blks
    
    # Initialize dictionary 
    ot_dab = {}
    
    # Shift the start of the trajectory to the begin variable 
    universe.trajectory[begin]
    
    # Keys for the dictionary above 
    sf_lbl = ["Avg Radius of gyration","Avg Sq. radius of gyration", "Avg end to end distance","Avg Sq. end to end distance", "Avg persistence length", "Avg Hydrodynamic radius"]

    # Array that will contain the Rg and Lp for each trajectory block 
    blk_nparr = np.zeros(shape=(len(sf_lbl),no_of_blks))
    
    # Arrays that will statistics of the Lp linear fitting 
    mod_res = np.zeros(shape=(4,no_of_blks))
    
    # initialize counter
    count = 0 

    for i in range(no_of_blks):
        
        # Shift the start of the trajectory to the start variable
        start = universe.trajectory.frame
        print(start)
    
        # Define end of block 
        end = int(start + n_size)
        print(end)
       
        # Get the time averaged cosine theta values for this block 
        rg_ens, cor_tp, theta_ens, rg_avg = get_rg_pers_poly(polymer_atoms, universe, start, end)
        
        # Using MD Analysis, I can get the number of polymer monomers in the polymer_atoms input variable. 
        n_monomers = len(np.unique(polymer_atoms.resids))
        
        # End to end distance and squared e2e distance
        eVec_poly, e2edis_poly = mean_sq_e2e(polymer_atoms, universe, n_monomers, start, end)
        
        # Store average radius of gyration for the block 
        blk_nparr[0,count] = rg_avg
        
        # Store averaged squared radius of gyration for the block 
        blk_nparr[1,count] = np.mean(rg_ens[1])
        
        # Store average end to end distance for the block 
        blk_nparr[2, count] = np.sqrt(np.mean(e2edis_poly[0]))
        
        # Store mean squared end to end distance for the block 
        blk_nparr[3, count] = np.mean(e2edis_poly[0])
        
        # Arc length x values
        blen = cor_tp[3]*len_bnd
        
        warnings.filterwarnings("error")

        try:
            # ln(cosine theta) y values
            np.log(cor_tp[0])
        except RuntimeWarning:
            
            print("Negative cosine theta values present") 
            
            # ln(cosine theta) y values
            npoly_lc = np.log(cor_tp[0])
        
            if fit_pnts == len(blen): 
            
                # Want to fit a line with no y-intercept 
                model_npoly = LinearRegression(fit_intercept=False)
            
                # fit line to data 
                model_npoly.fit(blen.reshape(-1,1), npoly_lc)

                # Predict new ln(cos(theta)) values from arc length values
                gg_np = model_npoly.predict(blen.reshape(-1,1))
            
                # Residuals between the true y data and model y data 
                resid_np = npoly_lc - gg_np
            
                # How to calculate mean squared error
                mse_p = np.sum(resid_np**2)/len(resid_np)
            
                # How to calculate Sum((Xi - avg(X))^2): X values are the arc length values 
                blen -= np.mean(blen)
                nhui = blen**2
                sum_m = np.sum(nhui)
            
                # How to calculate 95% confidence interval for the slope 
                flc_np = stats.t.ppf(0.975, fit_pnts-1)*np.sqrt(mse_p/sum_m)
            
                # Slope here is in angstroms
                print("Lp [Angstroms]:", -1/model_npoly.coef_[0])
            
                # Pers length error: error propagation from uncertainty in slope
                print("Error in Lp from fit [Angstroms], 95% CL:", flc_np/((model_npoly.coef_[0])**2))
            
                # Pearson coefficient to evaluate goodness of fit 
                print("R2 score:", sklearn.metrics.r2_score(npoly_lc, gg_np))
        
                # Save Lp in matrix
                mod_res[0,count] = -1/model_npoly.coef_[0]
            
                blk_nparr[4,count] = -1/model_npoly.coef_[0]
            
                # Save error in Lp from fit: Error propagation from the fit to Lp
                mod_res[1,count] = flc_np/((model_npoly.coef_[0])**2)
            
                # Save model slope 
                mod_res[2, count] = model_npoly.coef_[0]
            
                # Save Mean squared error of the fit 
                mod_res[3,count] = sklearn.metrics.mean_squared_error(npoly_lc, gg_np)
        
            elif fit_pnts != len(blen):
        
                # Want to fit a line with no y-intercept
                model_npoly = LinearRegression(fit_intercept=False)
            
                # fit line to data
                model_npoly.fit(blen[:fit_pnts].reshape(-1,1), npoly_lc[:fit_pnts])

                # Predict new ln(cos(theta)) values from arc length values
                gg_np = model_npoly.predict(blen[:fit_pnts].reshape(-1,1))
            
                # Residuals between the true y data and model y data 
                resid_np = npoly_lc[:fit_pnts] - gg_np[:fit_pnts]
            
                # How to calculate mean squared error
                mse_p = np.sum(resid_np**2)/len(resid_np)
            
                # How to calculate Sum((Xi - avg(X))^2): X values are the arc length values 
                blen -= np.mean(blen[:fit_pnts])
                nhui = blen**2
                sum_m = np.sum(nhui)
            
                # How to calculate 95% confidence interval for the slope 
                flc_np = scipy.stats.t.ppf(0.975, fit_pnts-1)*np.sqrt(mse_p/sum_m)
            
                # Slope here is in angstroms
                print("Lp [Angstroms]:", -1/model_npoly.coef_[0])
            
                # Pers length error: error propagation from uncertainty in slope
                print("Error in Lp from fit [Angstroms], 95% CL :", flc_np/((model_npoly.coef_[0])**2))
            
                # Pearson coefficient to evaluate goodness of fit
                print("R2 score:", sklearn.metrics.r2_score(npoly_lc[:fit_pnts], gg_np[:fit_pnts]))
        
                # Save Lp in matrix
                mod_res[0,count] = -1/model_npoly.coef_[0]
            
                blk_nparr[4,count] = -1/model_npoly.coef_[0]
            
                # Save error in Lp from fit: Error propagation from the fit to Lp
                mod_res[1,count] = flc_np/((model_npoly.coef_[0])**2)
            
                # Save model slope
                mod_res[2, count] = model_npoly.coef_[0]
            
                # Save Mean squared error of the fit 
                mod_res[3,count] = sklearn.metrics.mean_squared_error(npoly_lc[:fit_pnts], gg_np[:fit_pnts])
                
        else:
            
            # ln(cosine theta) y values
            npoly_lc = np.log(cor_tp[0])
        
            if fit_pnts == len(blen): 
            
                # Want to fit a line with no y-intercept 
                model_npoly = LinearRegression(fit_intercept=False)
            
                # fit line to data 
                model_npoly.fit(blen.reshape(-1,1), npoly_lc)

                # Predict new ln(cos(theta)) values from arc length values
                gg_np = model_npoly.predict(blen.reshape(-1,1))
            
                # Residuals between the true y data and model y data 
                resid_np = npoly_lc - gg_np
            
                # How to calculate mean squared error
                mse_p = np.sum(resid_np**2)/len(resid_np)
            
                # How to calculate Sum((Xi - avg(X))^2): X values are the arc length values 
                blen -= np.mean(blen)
                nhui = blen**2
                sum_m = np.sum(nhui)
            
                # How to calculate 95% confidence interval for the slope 
                flc_np = stats.t.ppf(0.975, fit_pnts-1)*np.sqrt(mse_p/sum_m)
            
                # Slope here is in angstroms
                print("Lp [Angstroms]:", -1/model_npoly.coef_[0])
            
                # Pers length error: error propagation from uncertainty in slope
                print("Error in Lp from fit [Angstroms], 95% CL:", flc_np/((model_npoly.coef_[0])**2))
            
                # Pearson coefficient to evaluate goodness of fit 
                print("R2 score:", sklearn.metrics.r2_score(npoly_lc, gg_np))
        
                # Save Lp in matrix
                mod_res[0,count] = -1/model_npoly.coef_[0]
            
                blk_nparr[4,count] = -1/model_npoly.coef_[0]
            
                # Save error in Lp from fit: Error propagation from the fit to Lp
                mod_res[1,count] = flc_np/((model_npoly.coef_[0])**2)
            
                # Save model slope 
                mod_res[2, count] = model_npoly.coef_[0]
            
                # Save Mean squared error of the fit 
                mod_res[3,count] = sklearn.metrics.mean_squared_error(npoly_lc, gg_np)
        
            elif fit_pnts != len(blen):
        
                # Want to fit a line with no y-intercept
                model_npoly = LinearRegression(fit_intercept=False)
            
                # fit line to data
                model_npoly.fit(blen[:fit_pnts].reshape(-1,1), npoly_lc[:fit_pnts])

                # Predict new ln(cos(theta)) values from arc length values
                gg_np = model_npoly.predict(blen[:fit_pnts].reshape(-1,1))
            
                # Residuals between the true y data and model y data 
                resid_np = npoly_lc[:fit_pnts] - gg_np[:fit_pnts]
            
                # How to calculate mean squared error
                mse_p = np.sum(resid_np**2)/len(resid_np)
            
                # How to calculate Sum((Xi - avg(X))^2): X values are the arc length values 
                blen -= np.mean(blen[:fit_pnts])
                nhui = blen**2
                sum_m = np.sum(nhui)
            
                # How to calculate 95% confidence interval for the slope 
                flc_np = scipy.stats.t.ppf(0.975, fit_pnts-1)*np.sqrt(mse_p/sum_m)
            
                # Slope here is in angstroms
                print("Lp [Angstroms]:", -1/model_npoly.coef_[0])
            
                # Pers length error: error propagation from uncertainty in slope
                print("Error in Lp from fit [Angstroms], 95% CL :", flc_np/((model_npoly.coef_[0])**2))
            
                # Pearson coefficient to evaluate goodness of fit
                print("R2 score:", sklearn.metrics.r2_score(npoly_lc[:fit_pnts], gg_np[:fit_pnts]))
        
                # Save Lp in matrix
                mod_res[0,count] = -1/model_npoly.coef_[0]
            
                blk_nparr[4,count] = -1/model_npoly.coef_[0]
            
                # Save error in Lp from fit: Error propagation from the fit to Lp
                mod_res[1,count] = flc_np/((model_npoly.coef_[0])**2)
            
                # Save model slope
                mod_res[2, count] = model_npoly.coef_[0]
            
                # Save Mean squared error of the fit 
                mod_res[3,count] = sklearn.metrics.mean_squared_error(npoly_lc[:fit_pnts], gg_np[:fit_pnts])
        
        count += 1
        
        # Shift the start to the next trajectory block 
        universe.trajectory[end]
        
    # Hydrodynamic Radius for this block
    rh_poly = rh_block_avg(no_of_blks, polymer_atoms, universe, begin, final)
    
    ot_dab[sf_lbl[5]] = rh_poly
   
    for i in range(len(sf_lbl)-1):
        ot_dab[sf_lbl[i]] = blk_nparr[i,:]
    
    return ot_dab, mod_res

# Name of paper: Computer simulation of dilute polymer solutions with the dissipative particle dynamics method 
# Authors: A. G. Schlijper, P. J. Hoogerbrugge, and C. W. Manke
def pos_bead_autocorr_RA(polymer_atoms, universe, n_monomers, t_corr, window_shift, start, end):
    
    """This function calculates the positional bead autocorrelation as a function of time, with running averaging method.
    t_corr is the number of frames in a window. """
    
    # Correlation vs time matrix
    pb_corr = np.zeros(shape=(int(t_corr)))
                       
    # Time lag array
    t_lag = np.arange(0, int(t_corr))
                       
    #Initialize matrix to store correlation values for each window
    # Default window displacement is one frame down
    #tcof_TO = np.zeros(shape=end-(t_corr-1), dtype=object)
                       
    #tcot_sum = np.zeros(shape=end-(t_corr-1), dtype=object)

    # counter for frame selection 
    #count = 0
    
    # Get No. of samples for a given window displacement
    n_wind = 0
    c_n = 0 
    e_n = start + t_corr

    for i in range(end):

        c_n += 1

        if c_n == e_n:
        
            n_wind += 1
        
            e_n = e_n + window_shift 

    print("No. of Samples: "+str(n_wind))

    #Initialize matrix to store correlation values for each window
    # Default window displacement is one frame down
    tcof_TO = np.zeros(shape=n_wind, dtype=object)
                       
    tcot_sum = np.zeros(shape=n_wind, dtype=object)

    for i in range(n_wind):
    
        if i == 0:
        
            # Shift the start of the trajectory to the start variable
            start_to = i
            #print(start_to)
        

            # Define end point based on no. of tie origins 
            end_to = int(t_corr)
            print(str(start_to)+ " to "+str(end_to))
        
            # Initialize matrix to store correlation values for each frame within each block 
            t_cofcorr = np.zeros(shape=end_to-start_to, dtype=object)
                       
            t_cofSUM = np.zeros(shape=end_to-start_to, dtype=object)
    
        elif i != 0:
        
            #print(i)
        
            # Shift the start of the trajectory to the start variable
            start_to = start_to + window_shift
            #print(start_to)
        
            # Define end point based on no. of tie origins 
            end_to = end_to + window_shift
            print(str(start_to)+ " to "+str(end_to))
        
            #Initialize matrix to store correlation values for each frame after the first block
            t_cofcorr = np.zeros(shape=end_to-start_to, dtype=object)
                       
            t_cofSUM = np.zeros(shape=end_to-start_to, dtype=object)
        
        # Initilaize variable for monomers COM
        com_OL = np.zeros(shape=(n_monomers), dtype=object)
    
        # New counter
        c_n2 = 0

        # Initialize var for olig com at the first time origin frame
        com_chain = []

        # For loop section, frame iteration 
        for ts in universe.trajectory[start_to:end_to]:
    
            # First frame(First time origin)
            if c_n2 == 0:
        
                # Initialize array for distance storage
                dist_com = np.zeros(shape=(n_monomers), dtype=object)
        
                # Initialize variable for dot product values
                n_bcorr = np.zeros(shape=(n_monomers))
        
                # Center of mass of the oligomer chain. save in list  
                com_chain.append(polymer_atoms.center_of_mass())
        
                for j in range(n_monomers):
            
                    #print(j+1)
            
                    a_fmon = polymer_atoms.select_atoms("resid "+str(j+1)).center_of_mass()
            
                    # Save COM for each monomer in an array
                    com_OL[j] = a_fmon
            
                    # save distances and normalize vector 
                    dist_com[j] = (a_fmon - com_chain[0])/(np.linalg.norm((a_fmon - com_chain[0])))
            
                    # Take dot product
                    n_bcorr[j] = np.dot(dist_com[j],dist_com[j])
    
                t_cofcorr[c_n2] = n_bcorr
                       
                t_cofSUM[c_n2] = np.sum(n_bcorr)/n_monomers
            
            # Following frames 
            elif c_n2 != 0:
        
                # Initialize array for distance storage
                dist_Afcom = np.zeros(shape=(n_monomers), dtype=object)
        
                # Initialize variable for dot product values
                n_Afcorr = np.zeros(shape=(n_monomers)) 
        
                # Center of mass of the oligomer chain 
                com_olig = polymer_atoms.center_of_mass()
        
                for j in range(n_monomers):
            
                    #print(j+1)
            
                    # COM for each monomer 
                    a_NFmon = polymer_atoms.select_atoms("resid "+str(j+1)).center_of_mass()
            
                    # save distances and normalize vector 
                    dist_Afcom[j] = (a_NFmon - com_olig)/(np.linalg.norm((a_NFmon - com_olig)))
            
                    # Take dot product
                    n_Afcorr[j] = np.dot(dist_Afcom[j],dist_com[j])
        
                t_cofcorr[c_n2] = n_Afcorr
                
                t_cofSUM[c_n2] = np.sum(n_Afcorr)/n_monomers
        
            #print(n6_plga_ace.trajectory.frame)
        
            c_n2 += 1
            #print(c_n2)
    
            #count += 1
            #print(count)
        
        # Save correlation data vs time for each block 
        tcof_TO[i] = t_cofcorr
        
        tcot_sum[i] = t_cofSUM
        
    
    #print(tcof_TO)
    
    #print(tcot_sum)
    
    # Initiallize array to store time averaged correlation values for each monomer, for each time lag point
    sv_MN = np.zeros(shape=(t_corr), dtype=object)

    # Time averaging dot product for each monomer then summing to get the correlation values of the polymer vs time
    # Iterate through each time lag, based on no. of time origins 
    for j in range(t_corr):
    
        # Initialize array to store time averaged corrletation values for each monomer, from each block 
        ct_mean = np.zeros(shape=n_monomers)
    
        #print("j ="+str(j))
    
        # Iterate through monomers
        for k in range(n_monomers):
        
            sv_mon = []
        
            #print("k ="+str(k))
            
            # Iterate through time origins 
            for i in range(tcof_TO.shape[0]):
                #print("i ="+str(i))
                
                #print(i)
                
                # Save each correlation values across time blocks, for each monomer 
                sv_mon.append(tcof_TO[i][j][k])
            
            #print(sv_mon)
            
            # Time averaging happens here 
            ct_mean[k] = np.mean(sv_mon)
    
        # Correlation values for each time lag is calculated here 
        pb_corr[j] = np.sum(ct_mean)/n_monomers
    
        # Save mean values, output #1 
        # row 1: Correlation values at t = 0, each column is the time-averaged correlation for each monomer              
        
        sv_MN[j] = ct_mean
    
    # Output #2
    corr_ot = np.array([pb_corr, t_lag])
    
    return corr_ot, tcot_sum

# Name of paper: Computer simulation of dilute polymer solutions with the dissipative particle dynamics method 
# Authors: A. G. Schlijper, P. J. Hoogerbrugge, and C. W. Manke 
def rouse_relax(x, tr_1, n_bonds):
    
    ## Change here for fitting 
    f = 6/(((n_bonds)**2) - 1)

    a1_sv = []
    cor_times = []
    
    for i in range(n_bonds):
        
        #print(i)
        
        atv_i = 4*((np.sin(((i+1)*np.pi)/(2*n_bonds)))**2)
        
        if i == 0:
            
            a1_sv.append(atv_i)
            
            trouse_i = (atv_i/atv_i)*tr_1
            
            #trouse_i = (asg_i/asg_i)*tr_1
            
            cor_times.append((1/atv_i)*np.exp(-x/trouse_i))
            
        elif i != 0:
            
            trouse_i = (a1_sv[0]/atv_i)*tr_1
            
            cor_times.append((1/atv_i)*np.exp(-x/trouse_i))
            
    return f*np.sum(cor_times)


def zimm_relax_fit(x, tr_1, h_s, n_bonds):
    
    ## Change here for fitting
    f = 6/(((n_bonds)**2) - 1)
    
    b_pl = 1 - (1.66*(h_s**0.78))
    
    sigma = -1.40*(h_s**0.78)
    
    cr_ls = []
    a1_sv = []
    cor_times = []
    
    for i in range(n_bonds):
        
        #print(i)
        
        atv_i = 4*((np.sin(((i+1)*np.pi)/(2*n_bonds)))**2)
        
        #print(atv_i)
        
        asg_i = atv_i*b_pl*(((i+1)/n_bonds)**sigma) 
        
        #print(asg_i)
        
        if i == 0:
            
            a1_sv.append(atv_i)
            
            cr_ls.append(asg_i)
            
            trouse_i = (asg_i/asg_i)*tr_1
            
            #cor_times.append((1/atv_i)*np.exp(-x/trouse_i))
            
            cor_times.append((1/asg_i)*np.exp(-x/trouse_i))
            
        elif i != 0:
            
            trouse_i = (cr_ls[0]/asg_i)*tr_1
            
            #cor_times.append((1/atv_i)*np.exp(-x/trouse_i))
            
            cor_times.append((1/asg_i)*np.exp(-x/trouse_i))
            
    return f*np.sum(cor_times)

def zimm_relax_func(x, tr_1, h_s, n_bonds):
    
    ## Change here for fitting
    f = 6/(((n_bonds)**2) - 1)
    
    b_pl = 1 - (1.66*(h_s**0.78))
    
    sigma = -1.40*(h_s**0.78)
    
    cr_ls = []
    a1_sv = []
    t_zimm = []
    cor_times = []
    
    for i in range(n_bonds):
        
        #print(i)
        
        atv_i = 4*((np.sin(((i+1)*np.pi)/(2*n_bonds)))**2)
        
        #print(atv_i)
        
        asg_i = atv_i*b_pl*(((i+1)/n_bonds)**sigma)
        
        #print(asg_i)
        
        if i == 0:
            
            a1_sv.append(atv_i)
            
            cr_ls.append(asg_i)
            
            trouse_i = (asg_i/asg_i)*tr_1
            
            t_zimm.append(trouse_i)
            
            #cor_times.append((1/atv_i)*np.exp(-x/trouse_i))
            
            cor_times.append((1/asg_i)*np.exp(-x/trouse_i))
            
        elif i != 0:
            
            trouse_i = (cr_ls[0]/asg_i)*tr_1
            
            t_zimm.append(trouse_i)
            
            cr_ls.append(asg_i)
            
            #cor_times.append((1/atv_i)*np.exp(-x/trouse_i))
            
            cor_times.append((1/asg_i)*np.exp(-x/trouse_i))
            
    return f*np.sum(cor_times), t_zimm, cr_ls