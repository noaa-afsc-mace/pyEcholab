# coding=utf-8

#     National Oceanic and Atmospheric Administration (NOAA)
#     Alaskan Fisheries Science Center (AFSC)
#     Resource Assessment and Conservation Engineering (RACE)
#     Midwater Assessment and Conservation Engineering (MACE)

#  THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN
#  AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE FURNISHED "AS IS."
#  THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS,
#  EMPLOYEES, AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
#  OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY
#  (1) FOR THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
#  SUPPORT TO USERS.
"""
.. module:: echolab2.simrad_signal_proc

    :synopsis: Functions for processing Simrad EK80 data. Based on work
               by Chris Bassett, Lars Nonboe Anderson, Gavin Macaulay,
               Dezhang Chu, and many others.


| Developed by:  Rick Towler   <rick.towler@noaa.gov>
| National Oceanic and Atmospheric Administration (NOAA)
| Alaska Fisheries Science Center (AFSC)
| Midwater Assesment and Conservation Engineering Group (MACE)
|
| Author:
|       Rick Towler   <rick.towler@noaa.gov>
| Maintained by:
|       Rick Towler   <rick.towler@noaa.gov>
"""

import numpy as np

def create_ek80_tx(raw_data, calibration, return_pc=False,
        return_indices=None, fast=False):
    '''create_ek80_tx returns an array representing the ideal
    EK80 transmit signal computed using the parameters in the raw_data
    and calibration objects.

    fast (bool) - Set to true to assume all pings share the same Tx
                  signal. The signal and effective pulse length will
                  be computed for the first ping and replicated for
                  all others.

    '''

#TODO: The fast keyword is useless at the moment because it is not implemented in
#      the call to ek80.get_calibration (and subsequent calls). Since the first
#      step is to get cal parms, this method is initially called without fast set
#      and tx signals and associated vars are created for every ping. Subsequent
#      calls will set fast, but since the tx data are already cached, it is
#      ignored. The point of the fast keyword was to limit the need to compute
#      the tx signal for every ping, saving some time and memory to hold the
#      presumably identical data.
#
#      In order to implement this correctly, the fast keyword has to be exposed
#      in ek80.get_calibration OR, better yet, it should be set automatically
#      if all the pings share the same relevant params which should be determined
#      when the bulk of the parameters are extracted in ek80_calibration.from_raw_data()
#      Testing will then need to ensure that the cached data are broadcast correctly
#      when computing results.

#    print("CREATE_EK80_TX fast:" + str(fast))

    # If we're not given specific indices, grab everything.
    if return_indices is None:
        return_indices = np.arange(raw_data.ping_time.shape[0])

    # Ensure that return_indices is a numpy array
    elif type(return_indices) is not np.ndarray:
        return_indices = np.array(return_indices, ndmin=1)

    #  create a dict to hold the calibration data parameters we need
    cal_parms = {'slope':None,
                 'pulse_duration':None,
                 'frequency_start':None,
                 'frequency_end':None,
                 'frequency':None,
                 'rx_sample_frequency':None,
                 'filters':None}

    # Check if we have already computed the Tx signal - This function is
    # called multiple times during the conversion process so we cache the
    # results in the calibration object to speed things up. The cached
    # data is normally discarded at the end of the conversion process to
    # ensure stale data is not used during another conversion.
    cal_has_tx_data = False
    if hasattr(calibration, '_tx_signal'):
        # If the user has not set the 'fast' keyword, we check if the existing tx data
        # has the same number of pings as the size of the return indices. This is not a
        # foolproof check but it is the best we can do. If fast is set, then we just have
        # to assume that echolab manages the cache correctly or if the user is controlling
        # the cache that they know what they are doing.
        if fast or (calibration._tau_eff.shape[0] == return_indices.shape[0]):
                cal_has_tx_data = True

    # Check if we have cached data
    if not cal_has_tx_data:
        # No - compute new Tx signals

        # Get the cal params we need
        for key in cal_parms:
            cal_parms[key] = calibration.get_parameter(raw_data, key, return_indices)

        # if this is CW data, we will not have the start/end params so we
        # set them to cal_parms['frequency']
        if cal_parms['frequency_start'] is None:
            cal_parms['frequency_start'] = cal_parms['frequency']
            cal_parms['frequency_end'] = cal_parms['frequency']

        # Iterate thru the return indices
        if fast:
            #  if the fast keyword is set, we assume all of the pings will
            #  be the same so we only compute the first Tx signal
            return_indices = np.array(return_indices[0:1])

        # create the return arrays
        n_return = return_indices.shape[0]
        tx_data = []
        tau_eff = np.empty(n_return, dtype=np.float32)
        rx_sample_frequency_dec = np.empty(n_return, dtype=np.float32)

        for idx, ping_index in enumerate(return_indices):

            #  get the theoretical tx signal
            t, y_t = ek80_chirp2(cal_parms['frequency_start'][idx],
                    cal_parms['frequency_end'][idx], cal_parms['slope'][idx],
                    cal_parms['pulse_duration'][idx],
                    cal_parms['rx_sample_frequency'][idx])

            #  apply the stage 1 and stage 2 filters
            y, rx_sample_frequency_decimated = filter_and_decimate(y_t,
                    cal_parms['filters'][idx], [1, 2],
                    cal_parms['rx_sample_frequency'][idx])

            #  compute effective pulse duration
            if return_pc and raw_data.pulse_form[idx] > 0:
                y_eff = np.convolve(y, np.flipud(np.conj(y))) / np.linalg.norm(y,2) ** 2
            else:
                y_eff = y
            ptxa = np.abs(y_eff) ** 2
            teff = np.sum(ptxa) / (np.max(ptxa) * rx_sample_frequency_decimated)

            #  store this ping's tx signal
            tx_data.append(y)
            tau_eff[idx] = teff
            rx_sample_frequency_dec[idx] = rx_sample_frequency_decimated

# TODO: Remove this replicating of data when vectorizing in pulse_compression method
#        if fast:
#            #  We're assuming all pings are the same. Replicate first Tx signal to all pings.
#            tx_data = [y] * n_return
#            tau_eff[1:] = teff

        # cache the results in the calibration object
        calibration._tx_signal = tx_data
        calibration._tau_eff = tau_eff
        #calibration._y_t = y_t
        calibration._rx_sample_frequency_decimated = rx_sample_frequency_dec
        calibration._pulse_duration = cal_parms['pulse_duration']
        
        return tx_data, tau_eff #, y_t
    else:
        #  return the cached data
        return calibration._tx_signal, calibration._tau_eff #, calibration._y_t


def compute_effective_pulse_duration(raw_data, calibration, return_pc=False,
        return_indices=None, fast=False):

    # Get the ideal transmit pulse which also returns effective pulse duration
    tx_data, tau_eff = create_ek80_tx(raw_data, calibration, return_pc=return_pc,
            return_indices=return_indices, fast=fast)

    return tau_eff


def filter_and_decimate(y, filters, stages, rx_sample_frequency):
    '''filter_and_decimate will apply one or more convolution and
    decimation operations on the provided data.

        y (complex) - The signal to be filtered
        filters (dict) - The filters dictionary associated with the signal
                         being filtered.
        stages (int, list) - A scalar integer or list of integers containing
                             the filter stage(s) to apply. Stages are applied
                             in the order they are added to the list.

    '''
    #  make sure stages is iterable - make it so if not.
    try:
        iter(stages)
    except Exception:
        stages = [stages]

    #  get the provided filter stages
    filter_stages = filters.keys()

    #  initialize rx_sample_frequency_decimated
    rx_sample_frequency_decimated = rx_sample_frequency

    # iterate through the stages and apply filters in order provided
    for stage in stages:
        #  make sure this filter exists
        if not stage in filter_stages:
            raise ValueError("Filter stage " + str(stage) + " is not in the " +
                    "supplied filters dictionary.")

        # get the filter coefficients and decimation factor
        coefficients = filters[stage]['coefficients']
        decimation_factor = filters[stage]['decimation_factor']

        # apply filter
        y = np.convolve(y, coefficients)

        # and decimate
        y = y[0::decimation_factor]

        #  compute the decimated sampling rate
        rx_sample_frequency_decimated /= decimation_factor

    return y, rx_sample_frequency_decimated


def calc_hanning_window(raw_data, calibration, return_indices=None, fast=False,
        **kwargs):

    '''calc_hanning_window computes the


        Length of Hanning window currently chosen as 2^k samples for
        lowest k where 2^k >= 2 * No of samples in pulse


       The function is largely copied from code provided with the following paper:

        Andersen, L. N., Chu, D. Heimvoll, H, Korneliussen, R, Macaulay, G, Ona, E.
        Patel R., & Pedersen G. (2021, Apr. 15). Quantitative processing of broadband data
        as implemented in a scientific splitbeam echosounder. ArXiv.
        https://doi.org/10.48550/arXiv.2104.07248


    fast (bool) - Set to true to assume all pings share the same Tx
                  signal. This keyword will be ignored if the provided
                  calibration object has cached tx signal attributes
                  since

    '''

    # If we're not given specific indices, grab everything.
    if return_indices is None:
        return_indices = np.arange(raw_data.ping_time.shape[0])

    # Ensure that return_indices is a numpy array
    elif type(return_indices) is not np.ndarray:
        return_indices = np.array(return_indices, ndmin=1)

    # Check if there is cached intermediate values and if they are the
    # correct dimensions
    cal_has_tx_data = False
    if hasattr(calibration, '_tx_signal'):
        if fast or (calibration._tau_eff.shape[0] == return_indices.shape[0]):
            cal_has_tx_data = True

    # No cached data, we need to compute it
    if not cal_has_tx_data:
        _, _ = create_ek80_tx(raw_data, calibration, return_indices=return_indices,
                fast=fast, **kwargs)

    #  specify the cal params we need and get them
    cal_parms = {'sample_interval':None,
                 'pulse_duration':None,
                 'sound_speed':None}
    for key in cal_parms:
        cal_parms[key] = calibration.get_parameter(raw_data, key, return_indices)

    # Iterate thru the return indices
    if fast:
        # If the fast keyword is set, we assume all of the pings will  be the same
        # so we only compute the hanning window for the first requested return index.

        #  set return_indices to the first requested index
        return_indices = np.array(return_indices[0], ndmin=1)

        # Compute the sample thickness
        thickness = cal_parms['sample_interval'][0] * cal_parms['sound_speed'][0] / 2.0

        # Compute the filter length using parameters from the first returned ping.
        L = (cal_parms['sound_speed'][0] * 2 * cal_parms['pulse_duration'][0]) / thickness

        # Compute the length of the window for each ping
        N_w = np.array(2 ** np.ceil(np.log2(L)).astype('int32'), ndmin=1)

    else:
        # Fast keyword is unset - we compute hanning windows and associated params for
        # all return indices.

        # Compute the sample thickness
        thickness = cal_parms['sample_interval'][0] * cal_parms['sound_speed'][0] / 2.0

        # Compute the filter length by ping.
        L = (cal_parms['sound_speed'] * 2 * cal_parms['pulse_duration'][0]) / thickness

        # Compute the length of the window for each ping
        N_w = 2 ** np.ceil(np.log2(L)).astype('int32')
        # or : N_w = np.ceil(2 ** np.log2(L))

    # create the return arrays
    n_return = return_indices.shape[0]
    w_tilde_i = []
    t_w = np.empty(n_return, dtype=np.float32)
    t_w_n = []

    for idx, ping_index in enumerate(return_indices):

        t_w_n.append(np.arange(0, N_w[idx]) / calibration._rx_sample_frequency_decimated[idx])
        t_w[idx] = N_w[idx] / calibration._rx_sample_frequency_decimated[idx]
        n_h = np.arange(0, N_w[idx])
        w_i =  0.5 * (1.0 - np.cos(2.0 * np.pi * n_h / (N_w[idx] - 1)))
        w_tilde_i.append(w_i / (np.linalg.norm(w_i) / np.sqrt(N_w[idx])))

    return w_tilde_i, N_w, t_w, t_w_n


def freqtransf(FFTvecin, fsdec, fvec):
        """
        Shift FFT frequencies for specified frequencies.
        
        Parameters
        ----------
        FFTvecin : np.array
            FFT data from decimated frequencies
        fsdec : float
            Decimated sampling frequency [Hz]
        fvec : np.array
            Specified frequencies [Hz]
                    
        Returns
        -------
        float
            Vector with corrected frequencies [Hz]
        """

        nfft = len(FFTvecin)
        idxtmp = np.floor(fvec / fsdec * nfft).astype("int")
        idx = np.mod(idxtmp, nfft)

        return FFTvecin[idx]

def calcPower(y_pc, z_td_e, z_rx_e, nu):
    """
    Calculate the received power into a matched load.
    
    Output received power values of 0.0 are set to 1e-20.
    
    Parameters
    ----------
    y_pc : np.array
        Pulse compressed signal [V]
    z_td_e : float
        Transducer electrical impedance [Ω]
    z_rx_e : float
        Receiver electrical impedance [Ω]
    nu : int
        Number of receiver channels [1]
    
    Returns
    -------
    np.array
        Received electrical power [W]
        
    """
    K1 = nu / ((2 * np.sqrt(2)) ** 2)
    K2 = (np.abs(z_rx_e + z_td_e) / z_rx_e) ** 2
    K3 = 1.0 / np.abs(z_td_e)
    C1Prx = K1 * K2 * K3
    Prx = C1Prx * np.abs(y_pc) ** 2
    Prx[Prx == 0] = 1e-20

    return Prx

def calcPulseCompSphericalSpread(y_pc_n, r_c_n):
    """
    Calculate the spherical spreading compensation.
    
    Parameters
    ----------
    y_pc_n : np.array
        Pulse compressed signal averaged over all transducer sectors [V]
    r_c_n : float
        Range to the centre of of the range volume covered by the sliding 
        window [m]
        
    Returns
    -------
    np.array
        Pulse compressed signal compensated for spherical spreading [Vm]
    """
    y_pc_s_n = y_pc_n * r_c_n

    return y_pc_s_n


def calc_DFT_for_Sv(calibration, y_pc_s_n, w_tilde_i, y_mf_auto_n, N_w,
                 r_c_n, step):

    # Prepare for append
    Y_pc_v_m_n = []
    Y_tilde_pc_v_m_n = []
    svf_range = []

    f_m = calibration.frequency_fft
    f_s_dec = calibration._rx_sample_frequency_decimated[0]

    # DFT of auto correlation function of the matched filter signal
    _Y_mf_auto_m = np.fft.fft(y_mf_auto_n, n=N_w)
    Y_mf_auto_m = freqtransf(_Y_mf_auto_m,
                                         f_s_dec, f_m)

    min_sample = 0  # int(r0 / dr)
    max_sample = len(y_pc_s_n)  # int(r1 / dr)

    bin_start_sample = min_sample
    bin_stop_sample = bin_start_sample + N_w
    n_bins = 0
    while (bin_stop_sample < max_sample):

        # Windowed data
        yspread_bin = w_tilde_i * y_pc_s_n[bin_start_sample:bin_stop_sample]

        # TODO: Consider calculating range values simply as (bin_stop_sample + bin_start_sample) / 2
        # Range for bin
        bin_center_sample = int((bin_stop_sample + bin_start_sample) / 2)
        bin_center_range = r_c_n[bin_center_sample]
        svf_range.append(bin_center_range)

        # DFT of windowed data
        _Y_pc_v_m = np.fft.fft(yspread_bin, n=N_w)
        Y_pc_v_m = freqtransf(_Y_pc_v_m, f_s_dec, f_m)

        # Normalized DFT of windowed data
        Y_tilde_pc_v_m = Y_pc_v_m / Y_mf_auto_m

        # Append data
        Y_pc_v_m_n.append([Y_pc_v_m])
        Y_tilde_pc_v_m_n.append([Y_tilde_pc_v_m])

        # Next range bin
        bin_start_sample += step
        bin_stop_sample = bin_start_sample + N_w
        n_bins += 1

    svf_range = np.array(svf_range)

    return Y_pc_v_m_n, Y_mf_auto_m, Y_tilde_pc_v_m_n, svf_range


def calcPowerFreqSv(Y_tilde_pc_v_m_n, N_u, z_rx_e, z_td_e):
        """
        Calculate the received power spectrum for the sliding window.
        
        Parameters
        ----------
        Y_tilde_pc_v_m_n : np.array
            DFT of the pulse compressed signal from a volume normalised by the DFT
            of the reduced autocorrelation function for the matched filter,
            compensated for spreading loss [1]
        N_u : int
            Number of transducer sectors/receiver channels
        z_td_e : float
            Transducer sector electric impedance [Ω]
        z_rx_e : float
            Receiver electric impedance [Ω]
        
        Returns
        -------
        np.array
            DFT of the received electric power in a matched load for the signal
            from a volume [Wm^2]
        """
        
        # Initialize list of power values by range
        P_rx_e_v_m_n = []

        # Impedances
        Z = (np.abs(z_rx_e + z_td_e) / np.abs(z_rx_e)) ** 2 / np.abs(z_td_e)

        # Loop over list of FFTs along range
        for Y_tilde_pc_v_m in Y_tilde_pc_v_m_n:
            P_rx_e_v_m = N_u * (np.abs(Y_tilde_pc_v_m) / (2 * np.sqrt(2))) ** 2 * Z
            # Append power to list
            P_rx_e_v_m_n.append(P_rx_e_v_m)

        return P_rx_e_v_m_n


def calcSvf(P_rx_e_t_m_n, alpha_m, p_tx_e, lambda_m, t_w,
                psi_m, g_0_m, c, svf_range):
        """
        Calculate Sv as a function of frequency.
        
        Parameters
        ----------
        P_rx_e_t_m_n : np.array
            DFT of the received electric power [W]
        alpha_m : float
            Acoustic absorption [dB/m]
        p_tx_e : float
            Transmitted electric power [W]
        lambda_m : float
            Acoustic wavelength [m]
        t_w : float
            Sliding window duration [s]
        psi_m : float
            Equivalent beam angle [sr]
        g_0_m : float
            Transducer gain [dB]
        c : float
            Speed of sound [m/s]
        svf_range : np.array
            Range [m]
            
        Returns
        -------
        np.array
            Sv(f) [dB re 1 m^-1]
        """
        
        # Initialize list of Svf by range
        Sv_m_n = np.empty([len(svf_range), len(alpha_m)], dtype=float)

        G = (p_tx_e * lambda_m**2 * c * t_w * psi_m * g_0_m**2) / (32 * np.pi**2)
        n = 0

        # Loop over list of power values along range
        for P_rx_e_t_m in P_rx_e_t_m_n:
            Sv_m = (
                10 * np.log10(P_rx_e_t_m)
                + 2 * alpha_m * svf_range[n]
                - 10 * np.log10(G)
            )

            # Add to array
            Sv_m_n[n, ] = Sv_m
            n += 1

        return Sv_m_n

def ek80_chirp2(f0, f1, slope, tau, fs):
    '''ek80_chirp2 returns a representation of the EK80 transmit signal
    as (time, amplitude) with a maximum amplitude of 1. This method was
    derived from code accompanying:

    Andersen, L. N., Chu, D. Heimvoll, H, Korneliussen, R, Macaulay, G, Ona, E.
    Patel R., & Pedersen G. (2021, Apr. 15). Quantitative processing of broadband data
    as implemented in a scientific splitbeam echosounder. ArXiv.
    https://doi.org/10.48550/arXiv.2104.07248

    f0 (float) - The starting Frequency of the tx signal obtained
                    from the raw file configuration datagram in Hz.
    f1 (float) - The starting Frequency of the tx signal obtained
                    from the raw file configuration datagram in Hz.
    slope (float) - The slope of the signal ramp up and ramp down
    tau (float) - The commanded transmit pulse length in s.
    fs (float) - the receiver A/D sampling frequency in Hz

    '''

    nsamples = int(np.floor(tau * fs))
    t = np.linspace(0, nsamples - 1, num=nsamples) * 1 / fs
    a = np.pi * (f1 - f0) / tau
    b = 2 * np.pi * f0
    y = np.cos(a * t * t + b * t)
    L = int(np.round(tau * fs * slope * 2.0))  # Length of hanning window
    w = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(0, L, 1) / (L - 1)))
    N = len(y)
    w1 = w[0:int(len(w) / 2)]
    w2 = w[int(len(w) / 2):-1]
    i0 = 0
    i1 = len(w1)
    i2 = N - len(w2)
    i3 = N
    y[i0:i1] = y[i0:i1] * w1
    y[i2:i3] = y[i2:i3] * w2

    #  normalize
    y[:] = y / np.max(y)

    return t, y


def ek80_chirp(txpower, fstart, fstop, slope, tau, z, rx_freq):
    '''ek80_chirp returns a representation of the EK80 transmit signal
    as (time, amplitude) with a maximum amplitude of 1.

    This code was originally written in MATLAB by
        Chris Bassett - cbassett@uw.edu
        University of Washington Applied Physics Lab

    txpower (float) - The transmit power in watts
    fstart (float) - The starting Frequency of the tx signal obtained
                    from the raw file configuration datagram in Hz.
    fstart (float) - The starting Frequency of the tx signal obtained
                    from the raw file configuration datagram in Hz.
    slope (float) - The slope of the signal ramp up and ramp down
    tau (float) - The commanded transmit pulse length in s.
    z (float) - The transducer impedance obtained from the raw file
                configuration datagram in ohms.
    rx_freq (float) - the receiver A/D sampling frequency in Hz

    '''

    # Create transmit signal
    sf = 1.0 / rx_freq
    a  = np.sqrt((txpower / 4.0) * (2.0 * z))
    t = np.arange(0, tau, sf)
    nt = t.shape[0]
    nwtx = int(2 * np.floor(slope * nt))
    wtx_tmp = np.hanning(nwtx)
    nwtxh = int(np.round(nwtx / 2))
    wtx = np.concatenate([wtx_tmp[0:nwtxh], np.ones((nt - nwtx)), wtx_tmp[nwtxh:]])
    beta = (fstop - fstart) * (tau ** -1)
    chirp = np.cos(2.0 * np.pi * (beta / 2.0 * (t ** 2) + fstart * t))
    y_tmp = a * chirp * wtx

    # The transmit signal must have a max amplitude of 1
    y = y_tmp / np.max(np.abs(y_tmp))

    return t, y


def calcAutoCorrelation(tx_signal):
    """
    Calculate the autocorrelation of the matched filter signal.
    
    Parameters
    ----------
    tx_signal : np.array
        Transmit signal
    Returns
    -------
    y_mf_auto_n : np.array
        Autocorrelation of the input filter [1]
    """
    
    y_mf_n_conj_rev = np.conj(tx_signal)[::-1]
    y_mf_twoNormSquared = np.linalg.norm(tx_signal, 2) ** 2
    y_mf_n_conj_rev = y_mf_n_conj_rev
    y_mf_twoNormSquared = y_mf_twoNormSquared

    # Calculate auto correlation function for matched filter
    tx_data_auto = np.convolve(tx_signal, y_mf_n_conj_rev) / y_mf_twoNormSquared

    return tx_data_auto


def pulse_compression(raw_data, calibration, return_indices=None, fast=False):
    '''pulse_compression applies a matched filter to the received signal using
    the simulated Tx signal as the filter template. This method will only apply
    the filter to FM pings. It does nothing to CW pings.


    fast (bool) - Set to True if calibration parameters for all FM pings are
                  the same. When True, the Tx signal generated for the first
                  FM ping will be used for all other FM ping. When False,
                  a unique Tx signal is generated for each FM ping.
    '''

    # If we're not given specific indices, grab everything.
    if return_indices is None:
        return_indices = np.arange(raw_data.ping_time.shape[0])
    # Ensure that return_indices is a numpy array
    elif type(return_indices) is not np.ndarray:
        return_indices = np.array(return_indices, ndmin=1)

    # Check if we have any fm pings to process
    is_fm = raw_data.pulse_form[return_indices] > 0
    if np.any(is_fm):
        # we do have fm pings - get the indices for those fm pings
        fm_pings = return_indices[is_fm]

        # get the simulated tx signal for these pings
        tx_signal, _ = create_ek80_tx(raw_data, calibration,
                return_indices=fm_pings, return_pc = True, fast=fast)

        # create a copy of the data to operate on and return
        p_data = raw_data.complex.copy()

#TODO: VECTORIZE THIS

        # apply match filter to these pings
        for p_idx, tx in enumerate(tx_signal):
            # create matched filter using tx signal for this ping
            tx_mf = np.flipud(np.conj(tx))
            ltx = tx_mf.shape[0]
            tx_n = np.linalg.norm(tx, 2) ** 2

            # apply the filter to each of the transducer sectors
            for q_idx in range(p_data[p_idx,:,:].shape[1]):
                # apply the filter
                filtered = np.convolve(p_data[p_idx,:,q_idx], tx_mf) / tx_n
                # remove filter delay
                p_data[p_idx,:,q_idx] = filtered[ltx-1:]


#        convolve = np.vectorize(np.convolve, signature='(n),(m)->(k)')
#convolve(np.eye(4), [1, 2, 1])
#array([[1., 2., 1., 0., 0., 0.],
#       [0., 1., 2., 1., 0., 0.],
#       [0., 0., 1., 2., 1., 0.],
#       [0., 0., 0., 1., 2., 1.]])

    else:
        # don't copy if we're not modifying the data
        p_data = raw_data.complex

    return p_data
