# write function to average across frames to give ensembled averaged cosine theta values at each N - 1 value
import MDAnalysis as mda
import statsmodels as stats
import math
import numpy as np 
import pandas
import sklearn

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
    
    # Add 1 to all values to begin sequence from 1, instead of 0 
    len_vec += 1

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
    
    # Initialize a row vector to store the radius of gyration at each frame 
    rg_ens = np.zeros(shape=(1,end-start))
    
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
        rg_ens[0,count_rg] = polymer_atoms.radius_of_gyration()
        
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
    avg_rg = np.mean(rg_ens)
    
    return  rg_ens, v_poly, corr_v, avg_rg 

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
    sf_lbl = ["Avg Radius of gyration", "Avg persistence length"]

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
        
        # Store average radius of gyration for the block 
        blk_nparr[0,count] = rg_avg
        
        # Arc length x values
        blen = cor_tp[3]*len_bnd
        
        # ln(cosine theta) y values
        npoly_lc = np.log(cor_tp[0])
        
        # Use sklearn to do fitting
        from sklearn.linear_model import LinearRegression
        
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
            flc_np = scipy.stats.t.ppf(0.975, fit_pnts-1)*np.sqrt(mse_p/sum_m)
            
            # Slope here is in angstroms
            print("Lp [Angstroms]:", -1/model_npoly.coef_[0])
            
            # Pers length error: error propagation from uncertainty in slope
            print("Error in Lp from fit [Angstroms], 95% CL:", flc_np/((model_npoly.coef_[0])**2))
            
            # Pearson coefficient to evaluate goodness of fit 
            print("R2 score:", sklearn.metrics.r2_score(npoly_lc, gg_np))
        
            # Save Lp in matrix
            mod_res[0,count] = -1/model_npoly.coef_[0]
            
            blk_nparr[1,count] = -1/model_npoly.coef_[0]
            
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
            
            blk_nparr[1,count] = -1/model_npoly.coef_[0]
            
            # Save error in Lp from fit: Error propagation from the fit to Lp
            mod_res[1,count] = flc_np/((model_npoly.coef_[0])**2)
            
            # Save model slope
            mod_res[2, count] = model_npoly.coef_[0]
            
            # Save Mean squared error of the fit 
            mod_res[3,count] = sklearn.metrics.mean_squared_error(npoly_lc[:fit_pnts], gg_np[:fit_pnts])        
        
        count += 1
        
        # Shift the start to the next trajectory block 
        universe.trajectory[end]
   
    for i in range(len(sf_lbl)):
        ot_dab[sf_lbl[i]] = blk_nparr[i,:]
    
    return ot_dab, mod_res