import numpy as np
from . import mask,line # used by apply_exclusion_mask()
from scipy import signal # used by smooth_Sv
import warnings

def noise_from_passive(P_data,thresh=20,min_range=5,run_mean_weights=np.array([.25,.5,.25])):
    """ takes a passive echolab object and computes noise in power [dB re 1W]
    for all the channels indicated.  This should only be used on passive data
    
    inputs:
        P_data - passive channel processed data object in power
        
        
        thresh - threshold to apply to remove outliers 
        
        (obs-median >thresh is exlcluded)
        
        min_range - ecxlude range less than this to avoid transmit pulse
        
        run_mean_weights = weights for running mean e.g. [.25,w .5,.25]
            The length of this should be an odd number
            if don't want a running mean shoudl be 1
    
    outputs:
    
        return echolab processed data object with noise_P to each entry in ch_dict
        with one value per ping in the object
    """
    # estimate noise for each ping
    if len(P_data.data.shape) ==3:
        data=P_data.data[:,:,P_data.range > min_range] 
        noise=np.apply_along_axis(get_noise,2,data,thresh)
    else:
        data=P_data.data[:,P_data.range > min_range] 
        noise=np.apply_along_axis(get_noise,1,data,thresh)

    P_data.add_object_attribute('noise_P',noise) # add noise as  newa data attribute
    if len(run_mean_weights) >1 : # if we want to smooth
        P_data.noise_P=running_mean(P_data.noise_P,run_mean_weights)
        
    return P_data

def running_mean(arr, weights=np.array([.25,.5,.25])):
    '''running mean
    computes running mean with edges padded with first/last value.
    this is what I'd think if as desired standard behavior
    
    inputs
    arr - 1-d nparray of interest
    array - np array of weights (normalized to 1 in the function)
    
    returns
    a running mean weighted array of same length as orignal array
    '''
    if len(arr.shape) > 1:
        new_arr = np.zeros((arr.shape[0], arr.shape[1] ))
        weights=weights/sum(weights) # normalize weights to sum to 1
        # Pad the beginning and the end of the array
        pad_width = len(weights) // 2
        for i in range(arr.shape[0]):
            padded_arr = np.pad(arr[i,:], (pad_width, pad_width), mode='edge')
            new_arr[i] = np.convolve(padded_arr, weights, mode='valid')
        # Calculate the running mean using convolution
        return new_arr
    else:
        weights=weights/sum(weights) # normalize weights to sum to 1
        # Pad the beginning and the end of the array
        pad_width = len(weights) // 2
        padded_arr = np.pad(arr, (pad_width, pad_width), mode='edge')  # 'edge' mode replicates the edge values
        # Calculate the running mean using convolution
        return np.convolve(padded_arr, weights, mode='valid')

def get_noise(x,thresh): 
    tmp =x[(x-np.median(x))<thresh] # trim to remove anything that might be noise
    return np_linearavg(tmp) # return noise computed in linear units and then then backtransformed

# define a helper function
def np_linearavg(x): # linear average of np array 
    return 10*np.log10(np.mean(10**(x/10)))

def Sv_noise_from_estimate(data_object):
    if data_object.data.ndim ==3:
        noise_hold_f = np.zeros(data_object.data.shape)
        noise_hold_power = data_object.to_power()
        for i in range(data_object.noise_power.shape[0]):
            tile_power = np.tile(data_object.noise_power[i],(data_object.data.shape[2],1)).transpose()
            noise_hold_power.data[i] = tile_power
        noise_hold = noise_hold_power.to_Sv()
        data_object.add_object_attribute('noise_Sv',noise_hold.data)
    else:
        tile_power = np.tile(data_object.noise_power,(data_object.data.shape[1],1)).transpose()
        noise_hold = data_object.to_power()
        noise_hold.data = tile_power
        noise_hold = noise_hold.to_Sv()
        data_object.add_object_attribute('noise_Sv',noise_hold.data)
    del noise_hold

def noise_correct(data_object): # linear average of np array 
    ''' subtracts two 10*log10 Sv values and noise corrects
    cases where raw_Sv<=noise_Sv are set to -999
    cases where 
    inputs x, y  - 10*log10(value)
    outputs x-y (transfomred to linear units) then backstransformed difference in 10*log10 units
    Useage  z=np_linearminus(x,y)'''
    diff= (10**(data_object.data/10))-(10**(data_object.noise_Sv/10))
    with warnings.catch_warnings(): # ignore warning as I know that log10 of negative numnbers will result in warning
        warnings.simplefilter("ignore")
        Sv=10*np.log10(diff) # negative values will be Sv
    Sv[np.isnan(Sv)]=-999 # set all Sv values that are Nan's because noise correction is negative to Sv of -999
    data_object.data=Sv
    return data_object


def SNR_threshold(data_obj,SNR_thresh):
    ''' compute SNR following De Robertis and Higginbottom eq 9-10
    Takes an Sv processed data object smooths it using an evenly weighted
    2-d matrix of adjacent samples via convolution. Convolution is done in linear domain
    results back-transformed into the
    
    inputs
        noise_corr_Sv - processed data object of noise corrected Sv
        active_data_noise - noise object corresponding to noise_corr_Sv 
        SNR_thresh - SNR threshold in dB
        
    Outputs
        noise_corr_Sv data processed object with SNR<SNR_thresh replaced with -999
    '''    
    SNR=data_obj.data-data_obj.noise_Sv #compute SNR
    data_obj.data[SNR<SNR_thresh]=-999
    return data_obj