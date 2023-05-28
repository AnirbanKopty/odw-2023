from pycbc.waveform import get_fd_waveform, get_td_waveform
from pycbc.frame import read_frame
from pycbc.filter import highpass, resample_to_delta_t, matched_filter
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.vetoes import power_chisq
from pycbc.events.ranking import newsnr
import pycbc
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


def index_peak_gen(nsnr_x, channel_name_arr, mass_arr, factor_method=True, factor=7, cutoff=4):
    """
    This function takes in `nsnr_x` `dict` which has a structure of `nsnr_x[cname][mass]`, finds peak indices.  

    The triggers can be taken in two ways based on the `cutoff` over which it will be considered as a potential trigger,
        - `factor_method=True`: it takes `mean+std*factor` to be cutoff
        - `factor_method=False`: it takes a `cutoff`

    Returns a `dict` containing the peak indices in the same structure as `nsnr_x`
    """
    # We consider peaks to be signal if it's above factor*sigma level (mean+std*factor)
    peak_idx = {}
    for cname in channel_name_arr:
        peak_idx[cname] = {}
        for mass in mass_arr:

            # #
            # if mass != 50:
            #     continue
            # #

            nsnr_x[cname][mass] = abs(nsnr_x[cname][mass])
            if factor_method:
                cutoff = np.mean(nsnr_x[cname][mass]) + np.std(nsnr_x[cname][mass])*factor
            else:
                pass        # taking the given cutoff value

            bool_arr = nsnr_x[cname][mass] > cutoff
            nsnr_x_peaks = nsnr_x[cname][mass][bool_arr]

            peak_idx[cname][mass] = []
            if isinstance(nsnr_x[cname][mass], pycbc.types.timeseries.TimeSeries):
                for peaks in nsnr_x_peaks:
                    peak_idx[cname][mass].append(np.where(nsnr_x[cname][mass].data == peaks)[0][0])
            else:
                for peaks in nsnr_x_peaks:
                    peak_idx[cname][mass].append(np.where(nsnr_x[cname][mass] == peaks)[0][0])

            # #? to eradicate nearby points and take end points (indices)
            peak_idxx = []
            peak_idxx.append(peak_idx[cname][mass][0])
            for i in range(len(peak_idx[cname][mass])-1):
                # print(np.diff(peak_idx[cname][mass]))
                if np.diff(peak_idx[cname][mass])[i] > 1:
                    peak_idxx.append(peak_idx[cname][mass][i])
                    peak_idxx.append(peak_idx[cname][mass][i+1])
            peak_idxx.append(peak_idx[cname][mass][-1])
            peak_idx[cname][mass] = peak_idxx
            
            # #? to take the indices of the maxima/peaks
            peak_idxx = []
            for i in range(0, len(peak_idx[cname][mass]), 2):
                try:
                    tmp_peak_idx = np.argmax(nsnr_x[cname][mass][peak_idx[cname][mass][i]:peak_idx[cname][mass][i+1]])
                except:
                    tmp_peak_idx = 0
                peak_idxx.append(peak_idx[cname][mass][i]+tmp_peak_idx)
            peak_idx[cname][mass] = peak_idxx
    
    return peak_idx

def plot_triggers(nsnr_x, channel_name_arr, mass_arr, peak_idx, snr_x, tof, factor, full_range=False, mass_val=None):
    #! This script can work only when will show correct windows and everything only if channels are not more than 2
    """
    It plots the `nsnr_x` and shows each peaks in zoomed view.
    
    Since the kernel was crashing, I added `full_range` option, if `True` plots for all masses, but if `False`, plots for only `mass_val`
    """

    for mass in mass_arr:
        if not(full_range):
            if mass!=mass_val:
                continue

        print(f"\n\nMass: {mass}")

        #? to check which channel has highest peak counts, it is needed to generate no of subplots
        pk_cnt = {}
        for cname in channel_name_arr:
            pk_cnt[cname] = len(peak_idx[cname][mass])
        pk_cnt_max = np.max(list(pk_cnt.values()))
        cname_maxpk = list(pk_cnt.keys())[list(pk_cnt.values()).index( pk_cnt_max )]
        #// print(pk_cnt)
        
        # Plot the new SNR timeseries
        fig1, ax1 = plt.subplots(figsize=[14, 4])
        # creating multiple subplots by 4 in each row
        four = 4
        n = len(peak_idx[cname_maxpk][mass])
        m = int(n/four)
        if n%four!=0:
            m += 1
        fig, ax = plt.subplots(m, four, figsize=[14,4*m], sharey=True)
        ax = ax.reshape(ax.size)

        for cname, color in zip(channel_name_arr, ['tab:blue', 'tab:orange']):
            ax1.plot(snr_x[cname][mass].sample_times, nsnr_x[cname][mass], label=cname)
            ax1.axhline(np.mean(nsnr_x[cname][mass])+np.std(nsnr_x[cname][mass])*factor, linewidth=0.5, color=color, linestyle='--')
            for i in range(len(peak_idx[cname][mass])):
                ax1.axvline(snr_x[cname][mass].sample_times[peak_idx[cname][mass][i]], color=color, linestyle='--', linewidth=0.5)

                # shade the region around each Hanford peak that could have a peak in Livingston if from
                # an astrophysical source
                if cname == cname_maxpk:
                    ptime = snr_x[cname][mass].sample_times[peak_idx[cname][mass]][i]
                    ax1.axvspan(ptime - tof, ptime + tof, alpha=0.2, color=color)
        
            # for the subplots showing zoomed views
            for iax, i in zip(np.atleast_1d(ax), range(len(np.atleast_1d(ax)))):
                if i < n:
                    iax.plot(snr_x[cname][mass].sample_times, nsnr_x[cname][mass], '.-', label=cname)
                    iax.axhline(np.mean(nsnr_x[cname][mass])+np.std(nsnr_x[cname][mass])*factor, linewidth=0.5, color=color, linestyle='--')
                    try:
                        # print(cname, peak_idx[cname][mass][i], snr_x[cname][mass].sample_times[peak_idx[cname][mass][i]])
                        iax.axvline(snr_x[cname][mass].sample_times[peak_idx[cname][mass][i]], color=color, linestyle='--', linewidth=0.5)
                    except:
                        pass

                    if cname == cname_maxpk:
                        iax.set_xlim( snr_x[cname][mass].sample_times[peak_idx[cname][mass][i]] - tof - 0.001,
                                snr_x[cname][mass].sample_times[peak_idx[cname][mass][i]] + tof + 0.001  )
                        ptime = snr_x[cname][mass].sample_times[peak_idx[cname][mass]][i]
                        iax.axvspan(ptime - tof, ptime + tof, alpha=0.2, color=color)

        ax1.set_title(f'NewSNR Timeseries for mass {mass}')
        ax1.grid()
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Re-weighted Signal-to-noise')
        ax1.legend()
        plt.show()


        for cname in channel_name_arr:
            print(f"\nChannel: {cname}")
            for idx in peak_idx[cname][mass]:
                print("We found a signal at {}s with SNR {}".format(snr_x[cname][mass].sample_times[idx],
                                                                    nsnr_x[cname][mass][idx]))





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
