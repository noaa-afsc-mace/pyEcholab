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
    tmp =x[(x-np.nanmedian(x))<thresh] # trim to remove anything that might be noise
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


def _smooth_Sv_2d(sv_data_2d, sample_thickness, pings_smooth=5, m_smooth=2.0):
    '''
    Helper function to smooth a 2D Sv array (pings x samples) in linear domain
    using normalized convolution to handle NaNs and edges without dilation.
    '''
    samples_smooth = max(1, int(round(m_smooth / sample_thickness)))
    weights = np.ones((pings_smooth, samples_smooth), dtype=np.float64)
    weights /= np.sum(weights)

    # Convert to linear domain (s_v)
    linear_data = 10.0 ** (sv_data_2d / 10.0)
    valid_mask = ~np.isnan(linear_data)
    linear_filled = np.where(valid_mask, linear_data, 0.0)

    # Normalized 2D convolution
    num = signal.convolve2d(linear_filled, weights, mode='same', boundary='symm')
    denom = signal.convolve2d(valid_mask.astype(np.float64), weights, mode='same', boundary='symm')

    with np.errstate(divide='ignore', invalid='ignore'):
        smoothed_linear = np.where(denom > 0, num / denom, np.nan)
        smoothed_sv = np.where(smoothed_linear > 0, 10.0 * np.log10(smoothed_linear), -999.0)
    smoothed_sv[denom == 0] = np.nan
    return smoothed_sv


def smooth_Sv(data_obj, pings_smooth=5, m_smooth=2.0):
    '''
    Smooths an Sv processed data object in the linear domain using normalized 2D convolution.
    For 2D data (CW), smooths across (pings, samples).
    For 3D data (FM), smooths each frequency slice in 2D on a frequency-by-frequency basis.
    '''
    sample_thickness = getattr(data_obj, 'sample_thickness', None)
    if sample_thickness is None or sample_thickness <= 0:
        if hasattr(data_obj, 'range') and len(data_obj.range) > 1:
            sample_thickness = abs(data_obj.range[1] - data_obj.range[0])
        else:
            sample_thickness = 1.0

    if data_obj.data.ndim == 2:
        data_obj.data = _smooth_Sv_2d(data_obj.data, sample_thickness,
                                      pings_smooth=pings_smooth, m_smooth=m_smooth)
    elif data_obj.data.ndim == 3:
        n_freq = data_obj.data.shape[0]
        indices = range(n_freq)
        for i in indices:
            data_obj.data[i] = _smooth_Sv_2d(data_obj.data[i], sample_thickness,
                                             pings_smooth=pings_smooth, m_smooth=m_smooth)
    return data_obj


def SNR_threshold(data_obj, SNR_thresh, exclude_val=-999, smooth=False,
                  pings_smooth=5, m_smooth=2.0, raw_data=None):
    '''
    Computes SNR following De Robertis and Higginbottom (2007) eqs 9-10.

    If smooth=True:
        Smooths Sv in linear domain using a 2D window (pings x meters) via normalized
        convolution to derive the SNR mask. For 3-D arrays (FM), convolution
        is performed in 2D on a frequency-by-frequency basis (each frequency slice is
        smoothed and masked independently across pings and range, never across frequencies).
    If smooth=False:
        Skips smoothing entirely and evaluates SNR point-by-point for maximum speed.
    '''
    # Fast path: skip smoothing entirely
    if not smooth:
        SNR = data_obj.data - data_obj.noise_Sv
        data_obj.data[SNR < SNR_thresh] = exclude_val
        return data_obj

    # If smooth...
    source_data = raw_data if raw_data is not None else data_obj.data

    sample_thickness = getattr(data_obj, 'sample_thickness', None)
    if sample_thickness is None or sample_thickness <= 0:
        if hasattr(data_obj, 'range') and len(data_obj.range) > 1:
            sample_thickness = abs(data_obj.range[1] - data_obj.range[0])
        else:
            sample_thickness = 1.0

    if source_data.ndim == 2:
        # 2D data (CW): smooth across (pings, samples)
        smoothed_sv = _smooth_Sv_2d(source_data, sample_thickness,
                                    pings_smooth=pings_smooth, m_smooth=m_smooth)
        snr = smoothed_sv - data_obj.noise_Sv
        data_obj.data[snr < SNR_thresh] = exclude_val

    elif source_data.ndim == 3:
        # 3D data (FM/broadband): shape is (n_frequencies, n_pings, n_samples)
        # Smooth and threshold on a frequency-by-frequency basis (smooth per 2D per slice)
        n_freq = source_data.shape[0]
        indices = range(n_freq)

        for i in indices:
            raw_slice = source_data[i, :, :]
            noise_slice = data_obj.noise_Sv[i, :, :]

            smoothed_sv_i = _smooth_Sv_2d(raw_slice, sample_thickness,
                                          pings_smooth=pings_smooth, m_smooth=m_smooth)
            snr_i = smoothed_sv_i - noise_slice
            mask_i = (snr_i < SNR_thresh)
            data_obj.data[i, mask_i] = exclude_val

    else:
        SNR = data_obj.data - data_obj.noise_Sv
        data_obj.data[SNR < SNR_thresh] = exclude_val

    return data_obj