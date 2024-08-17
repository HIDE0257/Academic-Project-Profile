from gc import collect
from math import radians, sin, cos, tan, floor, ceil, atan2, sqrt, pi, exp
from umatrix import zeros, matrix  # Matrix handler for micropython (refer to https://github.com/iyassou/umatrix)
# Multithreading
from _thread import start_new_thread

# Logging fn
from log import log_and_print

# Check out gnssrefl.gps.py.strip_compute in Dr. Larson's code

def reflection_height(wl,time,e_deg,az_deg,snr_dB,minH,maxH,min_az,max_az,min_e,max_e, corr_bool, press, temp_C, hum):    # Returns dt, H_bar, de

    e_rad = [radians(i) for i in e_deg]    # [rad] Elevation angle
    # Free the memory from e_deg
    del e_deg

    # Only run the computation if the range of the elevation angles is more than 10 degrees
    if(max(e_rad) - min(e_rad) < 0.174):
        #log_and_print("Range of Elevation angles is <10 degrees for this curve, exiting")
        return 0, 0, 0

    n = len(time)
    if n > 2 :
        # Flatten snr curve
        snr_linear = [10**(i/20) for i in snr_dB]
        snr_flat = flat_snr(snr_linear, time)   # [dB] flattened SNR corresponding to time
        if snr_flat == 0:
            return 0, 0, 0
        # Collect garbage from this fn
        del snr_dB, snr_linear
        collect()

        # Fit curve to elevation angle
        e_cont = fit_e(e_rad, time)
        # Collect garbage from this fn
        collect()
        
        # Mask the data here 
        idx_mask = masking(e_cont,az_deg,min_az,max_az,min_e,max_e)
        del az_deg, min_az, max_az, min_e, max_e
        collect()

    else:
        #log_and_print('The sampling time is not enough for this curve, exiting')
        return 0, 0, 0
    
    
    if len(idx_mask) > 35:
        # Below here, all variables are masked and are in a list
        time = [time[i] for i in idx_mask]
        dt = (time[-1] - time[0])/2
        snr_VV = [snr_flat[i] for i in idx_mask]    # This SNR is in V/V and flaten
        del time, snr_flat
        e_rad = [e_cont[i] for i in idx_mask]

            # Only run the computation if the range of the elevation angles is more than 10 degrees
        if(max(e_rad) - min(e_rad) < 0.174):
            #log_and_print("Range of Elevation angles is <10 degrees for this curve, exiting")
            return 0, 0, 0

        # Apply elevation angle correction, if we have that information
        if corr_bool:
            temp_K = temp_C + 273.15
            # Apply dynamic height correction
            # Elevation angle correction constants
            # Calculating wet pressure, P_w (refer to https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml)
            sat_press = 6.1094 * exp((17.625*temp_C)/(temp_C + 243.04)) # [hPa]
            P_w = sat_press*(hum/100)   # [hPa]
            N_0 = (77.689 * press / temp_K) + (71.2952 * P_w/temp_K) + (375463 * P_w/(temp_K**2))
            # Correct and replace the elevation angle
            for i in range(len(e_rad)):
                e_rad[i] = e_rad[i] + (10**-6) * N_0* cos(e_rad[i]) / (sin(e_rad[i]) + (0.00175*tan(1.527163095-e_rad[i]))) # [rad]
            del sat_press, P_w, N_0, temp_C, temp_K, press, hum
        else:
            pass
        collect()

        t = [sin(i)*(2/wl) for i in e_rad]
        de = (max(e_rad) - min(e_rad))/2
        # Garbage collection
        del e_rad, idx_mask, wl
        collect()

        # Get oversampling and high-frequency fators 
        ofac, hifac, RH = get_ofac_hifac_RH(t, maxH, 0.025)

        if ofac == 0 or hifac == 0 or RH == 0:
            return 0,0,0
        del ofac, hifac
        collect()

        # Now move onto LSP --< returns max reflector height [m]
        H_bar = LSP(t,snr_VV,RH,minH)
        #log_and_print(f"RH is {round(H_bar,2)} m,   dt = {round(dt,2)} s,    de = {round(de,2)} rad\n")

    else :
        #log_and_print('There are no datapoints after masking')
        return 0,0,0

    return int(dt), H_bar, de


def flat_snr(snr, time) : # Returns flaten snr
    # Inputs:   snr --> n-vector SNR [dB]
    #           time --> n-vector of time [s]
    # Outputs:  flat_snr --> n-vector of flatten snr [dB]
    #           e_hat --> [a;b;c] for e vs time 

    n = len(time)  # Length of array

    # We want quadratic fit (e(t) = a*t^2 + b*t + c --> e(t) = A*e_hat)
    #                       (snr(t) = a*t^2 + b*t + c --> snr(t) = A*x_hat)
    # Initialize an A matrix A = [t^2 t 1] 
    A = zeros(n,3)

    # Build an A matrix (n,3)
    for i in range(n):
        A[i][0:3] = [time[i]**2, time[i], 1]

    # Solve for x_hat --> [a; b; c]
    snr_mat = matrix(snr).transpose # SNR (n,1) in a matrix form
    try:
        x_hat = (A.transpose * A).inverse * (A.transpose*snr_mat)
    except Exception as e:
        log_and_print(e)
        return 0

    # Predict snr (Check if x_hat is a col vec or a row vec)

    if x_hat.shape[0] == 3 :
        snr_pred = A*x_hat
    else:
        snr_pred = A*x_hat.transpose
    del A, x_hat
    collect()

    # Take a residual to flat the SNR
    snr_flat = zeros(n,1)
    for i in range(n):
        snr_flat[i][0] = snr_mat[i][0] - snr_pred[i][0]
    
    # Convert matrix to list
    snr_flat = [i[0] for i in snr_flat]

    return snr_flat    # List respectively

def fit_e(e_rad, time) : # Returns flaten snr
    # Inputs:   e --> n-vector elevaton [rad]
    #           time --> n-vector of time [s]
    # Outputs:  flat_snr --> n-vector of flatten snr [dB]
    #           e_hat --> [a;b;c] for e vs time 

    n = len(time)  # Length of array

    # We want quadratic fit (e(t) = a*t^2 + b*t + c --> e(t) = A*e_hat)
    #                       (snr(t) = a*t^2 + b*t + c --> snr(t) = A*x_hat)
    # Initialize an A matrix A = [t^2 t 1] 
    A = zeros(n,3)

    # Build an A matrix (n,3)
    for i in range(n):
        A[i][0:3] = [time[i]**2, time[i], 1]

    # Solve for x_hat --> [a; b; c]
    e_mat = matrix(e_rad).transpose # SNR (n,1) in a matrix form
    try:
        e_hat = (A.transpose * A).inverse * (A.transpose*e_mat)
    except Exception as e:
        log_and_print(e)
        return 0
    collect()

    # Predict snr (Check if x_hat is a col vec or a row vec)
    if e_hat.shape[0] == 3 :
        e_fit = A*e_hat
    else:
        e_fit = A*e_hat.transpose
    del A, e_hat

    # Convert matrix to list
    e_fit = [i[0] for i in e_fit]

    return e_fit    # List respectively

def masking(e_cont,az_deg,min_az,max_az,min_e,max_e):
    n = len(e_cont)   # Length of array
    max_e_rad = radians(max_e)
    min_e_rad = radians(min_e)

    idx_mask = [None] * n
    for i in range(n):
        if min_az <= az_deg[i] and az_deg[i] <= max_az and min_e_rad <= e_cont[i] and e_cont[i] <= max_e_rad:
            idx_mask[i] = int(i)
    idx_mask = [i for i in idx_mask if i is not None]
        
    return idx_mask # List

def get_ofac_hifac_RH(t,maxH,precision) :
    # Inputs:   t = sin(e)*2/wl
    #           maxH = maximum reflector height (reflector height + 3m) [m] 
    #           precision = desired precision of Periodogram width [m]
    # Outputs:  ofac = oversampling factor 
    #           hifac = high-frequency factor
    #           RH = Reflector height array [m]
    #           These are used to define the range of frequency when computing LSP
   
    n = len(t)  # Number of observations
    dx = max(t) - min(t)   # Window length - difference between max and min of t

    if dx == 0:
        #log_and_print('Window length (x) is 0 - There is no datapoints after masking')
        return 0,0,0,0
    
    ofac = 1/(dx*precision)    # Oversampling factor 
    fc = n/(2*dx)   # [m] Nyquest Frequency  
    hifac = maxH/fc # high-freq factor 

    # Spacing out RH, using ofac and hifac
    rh0 = 1/(dx*ofac)   # [m] Reflector height spacing at start 
    rhn = hifac*n/(2*dx)# [m] Reflector height spacing at end
    m = int(0.5*ofac*hifac*n)   # Number of steps
    if m <= 1:
        return ofac,hifac,0,0
    else:
        steps = (rhn - rh0)/(m-1)   # [m] Increment 
        #log_and_print(f"m = {m}, n = {n}")

    RH = [None] * m     # [m] Initialize reflector height array
    # Make a linspace of fs
    for i in range(m):
        RH[i] = rh0 + i*steps # [m] Spacing out RH
    RH = [i for i in RH if i is not None]

        
    return ofac, hifac, RH  # Float, Float, List

# Implementation of LSP nested loop with multithreading
def LSP_nested(t, snr, m, r, RH, power_scaled):
    n = len(snr)

    for i in range(r[0], r[1]): 
        # Define the variables before running a loop
        tau_num = 0
        tau_den = 0
        pwr1_num = 0
        pwr1_den = 0
        pwr2_num = 0
        pwr2_den = 0

        w = 2*pi*RH[i] # Angular frequency

        for j in range(n):
            tau_num += sin(2*w*t[j])   # Tau numerator
            tau_den += cos(2*w*t[j])   # Tau denomerator

        # Time offset
        if tau_den == 0:
            tau = 0
        else:
            tau = (1/(2*w))*atan2(tau_num, tau_den)

        for k in range(n) :
            pwr1_num += snr[k]*cos(w*(t[k] - tau))   # pwr1 numerator
            pwr1_den += (cos(w*(t[k] - tau)))**2     # pwr1 denomerator

            pwr2_num += snr[k]*sin(w*(t[k] - tau))    # pwr2 numerator
            pwr2_den += (sin(w*(t[k] - tau)))**2     # pwr2 denomerator

        # Combine them above 
        pwr1 = (pwr1_num**2)/pwr1_den
        pwr2 = (pwr2_num**2)/pwr2_den
        # Define power [unitless]
        power = 0.5*(pwr1 + pwr2)

        # Scale the power 
        power_scaled[i + r[2]] = 2*sqrt(power/m)

def LSP(t,snr,RH,minH) :
    # Inputs:   t = sampled time array [s]
    #           snr = flattened snr array [V/V]
    #           RH = spaced-out reflection height array [m]
    #           minH = minimum height that we set [m]
    # Outputs:  maxP = maximum amplitude that gets from the LSP calculation
    #           maxRH = Reflector Height corresponding to the maxP [m]

    n = len(snr)    # Length of sample data
    m = len(RH)     # Span of RH

    # Predefine size of result array
    power_scaled = [1] * m

    # Threaded call to nested loop
    # Attempt to multithread, but fall back to core 0 if it fails
    bool_thread = True
    try:
        start_new_thread(LSP_nested, (t, snr, m, (0, ceil(m/2), floor(m/2)), RH[floor(m/2):], power_scaled))
    except OSError:
        # Thread failed, fallback to do the routine on core0
        bool_thread = False
    # process that is always on core 0
    LSP_nested(t, snr, m, (0, floor(m/2), 0), RH[:floor(m/2)], power_scaled)
    collect()
    if not bool_thread:
        LSP_nested(t, snr, m, (0, ceil(m/2), floor(m/2)), RH[floor(m/2):], power_scaled)

    # Collect garbage
    collect()

    # Trim anything less than minH
    for i in range(m):
        if RH[i] <= minH:
            power_scaled[i] = None
            RH[i] = None
    # Remove the none elements from the lists
    power_scaled = [i for i in power_scaled if i is not None]
    RH = [i for i in RH if i is not None]

    if len(power_scaled) == 0:
        #log_and_print('Invalid LSP, no data, check your inputs, make sure the masks are correct')
        maxRH = 0
        maxP = 0
    else:
        maxP = max(power_scaled)    # Maximum amplitude (power) Jan was here:) :/
        idx = argmax(power_scaled,maxP) # returns index of maxP
        if idx == -1:
            return 0
        else:
            maxRH = RH[idx]    # [m] Reflection height corresponding to maxP

    return maxRH   # Float

def argmax(array,lookfor):  # To find an index that corresponds to a x or y value
    # Inputs:   array = an array 
    #           lookfor = a value that you look for
    # Output:   idx_lookfor = index of the value

    for i in range(len(array)):
        if lookfor == array[i]:
            return i  # Integer

    #log_and_print("The argmax() doesn't work, check your trimed power in the LSP function")
    return -1
