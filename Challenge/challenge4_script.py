from pycbc.waveform import get_fd_waveform, get_td_waveform
from pycbc.frame import read_frame
from pycbc.filter import highpass, resample_to_delta_t, matched_filter
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.vetoes import power_chisq
from pycbc.events.ranking import newsnr
import matplotlib.pyplot as plt
import numpy as np

def challenge_matched_filter(file_name, channel_name_arr, mass_arr, spin=0):
    """
    Generalized to take single or multiple channels, and masses  

    Returns `ts`, `psd`, `snr`, `chisq`, `nsnr` `dict` corresponding to `channel_name` and `mass`
    """
    # print("\nLooking at file {} with template mass {} M_sol and spin {}".format(file_name, mass, spin))
    channel_name_arr = np.atleast_1d(channel_name_arr)
    mass_arr = np.atleast_1d(mass_arr)
    ts = {}
    psd_ts = {}
    snr_x = {}
    chisq_x = {}
    nsnr_x = {}
    # peak_idx = {}
    for cname in channel_name_arr:
        ts[cname] = read_frame(file_name, cname)
        ts[cname] = highpass(ts[cname], 15.0)
        ts[cname] = resample_to_delta_t(ts[cname], 1.0/1024)    #for saving time
        ts[cname] = ts[cname].crop(2, 2)
        psd_ts[cname] = ts[cname].psd(4)
        psd_ts[cname] = interpolate(psd_ts[cname], ts[cname].delta_f)
        psd_ts[cname] = inverse_spectrum_truncation(psd_ts[cname], int(4 * ts[cname].sample_rate),
        low_frequency_cutoff=15)
        
        snr_x[cname] = {}
        chisq_x[cname] = {}
        nsnr_x[cname] = {}
        # peak_idx[cname] = {}
        for mass in mass_arr:
            hp_x, _ = get_fd_waveform(approximant="IMRPhenomD",
                                    mass1=mass, mass2=mass, spin=spin,
                                    f_lower=20.0, delta_f=ts[cname].delta_f)
            hp_x.resize(len(psd_ts[cname]))

            # For each observatory use this template to calculate the SNR time series
            snr_x[cname][mass] = matched_filter(hp_x, ts[cname], psd=psd_ts[cname], low_frequency_cutoff=20).crop(5, 4)

            nbins = 26
            chisq_x[cname][mass] = power_chisq(hp_x, ts[cname], nbins, psd_ts[cname], low_frequency_cutoff=20.0)
            chisq_x[cname][mass] = chisq_x[cname][mass].crop(5, 4)

            dof_x = nbins * 2 - 2
            chisq_x[cname][mass] /= dof_x


            # The rho-hat term above is named "newsnr" here
            nsnr_x[cname][mass] = newsnr(abs(snr_x[cname][mass]), chisq_x[cname][mass])
    
    return ts, psd_ts, snr_x, chisq_x, nsnr_x


# Now that we've calculated the onsource peak, we should calculate the background peak values.
# We do this by chopping up the time series into chunks that are the same size as our
# on-source window and repeating the same peak finding (max) procedure - keeping the algorithm
# the same to prevent bias

# Walk through the data in chunks and calculate the peak statistic value in each.
def p_value_calculator(nsnr, sidx, window_size, plot=False, channel_name=None, mass=None):
    """
    It calculates p_value by chopping up the time series into chunks that are the same size as our 
    on-source window and repeating the same peak finding (max) procedure - keeping the algorithm
    the same to prevent bias

    ## Parameters
    `nsnr` : re-weighted statsitical snr value (which denies some glitches)  

    `sidx` : starting indeices of the peaks (can be an array for multiple peaks)  

    `window_size` : the size of the window around the peaks, use time of flight window, note that 
    this is in units of indices, not in time i.e., multiply by `sample_rate`, and it needs to be `int`  

    `plot` (bool) : if plotting visualization is needed, by default it is False  

    `channel_name` : when plot=True, this will help to identify the plots (plot title)  

    `mass` : same as `channel_name`

    Returns `dict` for p_values corresponding to `channel_name` and `mass`
    """
    # Calculate the span of time that a Virgo peak could in principle happen in from time of flight
    # considerations.
    # start = ptime['H1'] - tof['H1']
    # end = ptime['L1'] + tof['L1']

    # convert the times to indices along with how large the region is in number of samples
    sidx = np.atleast_1d(sidx)
    eidx = sidx + window_size

    # Calculate the "on-source" peak re-weighted (newsnr) statistic value.
    onsource = []
    for i in range(len(sidx)):
        onsource.append(nsnr[sidx[i]:eidx[i]].max())

    peaks = []
    i = 0
    while i + window_size < len(nsnr):
        p = nsnr[i:i+window_size].max()
        peaks.append(p)
        i += window_size
        
        # Skip past the onsource time
        for sidx_i in sidx:
            if abs(i - sidx_i) < window_size:
                i += window_size * 2
        
    peaks = np.array(peaks)


    # The p-value is just the number of samples observed in the background with a 
    # value equal or higher than the onsource divided by the number of samples.
    # We can make the mapping between statistic value and p-value using our background
    # samples.
    peaks.sort()

    pvalue = []
    for i in range(len(sidx)):
        pvalue.append((peaks > onsource[i]).sum() / float(len(peaks)))

    if plot:
        pcurve = np.arange(1, len(peaks)+1)[::-1] / float(len(peaks))
        plt.figure(figsize=[6, 5])
        plt.scatter(
            peaks, pcurve, label='Off-source (Noise Background)', color='black', s=1)

        for i in range(len(onsource)):
            plt.axvline(onsource[i], label='On-source', color='red')
            plt.axhline(pvalue[i], color='red')

        plt.legend()
        plt.yscale('log')
        plt.grid()
        plt.title(f"p-value estimation for channel {channel_name} with mass {mass}")
        # plt.ylim(1e-3, 1e0)
        plt.ylabel('p-value')
        plt.xlabel('Re-weighted Signal-to-noise')
        # plt.xlim(2, 5)
        plt.show()

    return pvalue
