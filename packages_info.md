# Detectors
- `'G1'` - GEO600
- `'H1'` - LIGO-Hanford
- `'L1'` - LIGO-Livingston
- `'V1'` - (Advanced) Virgo
- `'K1'` - KAGRA
- `'I1'` - Indigo?

# Packages for LIGO Data Analysis
- gwosc
- gwpy

## `gwosc` (0.7.1)
 
mostly used for quering available datasets
- `datasets`  
[
    *general notes:*   
        Event name 'GW150914-v3' means GW detection on  2015-09-14 and it's of version 3;  
        Run name 'O2_16KHZ_R1' means O2 run with 16KHZ sampling rate and this is Release 1;
        GPS time, the time system here uses GPS time system, counts number of seconds from the start of GPS epoch at 1980-01-06 00:00. GWOSC provides [GPS time converter](https://www.gw-openscience.org/gps/") or use `gwpy.time`.
]
    - `find_datasets`  
    tools for searching datasets including events, catalogs and full run strain data releases
        - `type` : `'catalog'`, `'events'`, `'run'`  
        for getting data based on types
        - `catalog`  
        to get data of some type from a catalog  
        [Event Portal](https://gw-openscience.org/eventapi)
        - `detector`  
        to get data of some type from a detector
        - `segment` (expects list/tuple of two entries - gps time in secs)
        to get data of some type for particular time segments
        - `version`  
        if there are multiple versions of the same event, version can be specified by this argument
        - `host`  
        URL of host to query datasets, defaults to `https://gwosc.org`
    - `event_gps`  
    for quering the GPS time of a specific event, expects a string of Event name (no need to specify the version no of the event name) 
    - `event_at_gps`  
    for quering the event which occured at the given gps time
    - `run_segment`
    for getting the time segment for a particular run
    - `run_at_gps`
    for getting run based on given gps time
    - `query_events`   
    for quering events based on physical parameter criterions
    `query_events(select=["25 <= network-matched-filter-snr <= 30", "total-mass-source >= 30"])`
        - `select`  
            [Parameters:](https://www.gwosc.org/apidocs/#event5)
            - ``gps-time``,
            - ``mass-1-source``,
            - ``mass-2-source``,
            - ``network-matched-filter-snr``,
            - ``luminosity-distance``,
            - ``chi-eff``,
            - ``total-mass-source``,
            - ``chirp-mass``,
            - ``chirp-mass-source``,
            - ``redshift``,
            - ``far``,
            - ``p-astro``,
            - ``final-mass-source``

- `locate`
    - `get_event_urls`  
    provides urls for a particular event, we can fetch even a particular url based on criterias
        - `event`  
        for specifying events
        - `duration`
        duration in seconds
        - `detector`  
        specifying detector (see later for detector specifications)


## `gwpy` (3.0.4)
This is a package based on OOP. Mainly for handling and analyzing the GW data.
[a programming/Python note: all functions use snake_case; all class objects use camelCase or PascalCase]

- `timeseries`
    - `TimeSeries`  
    TimeSeries object that has 1D data of GW signal for a particular time segment for a particular detector.  
        - `.fetch_open_data`  
        Downloads data directly from https://www.gwosc.org and creates a TimeSeries object.
        - `.read`  
        reads TimeSeries data from a file  
        Has several methods
        - `.plot`
        - `.times`
        - `.fft`
        - `.asd`
        - `.spectrogram`
        - `.spectrogram2`
            - `fftlength`
            - `overlap`
            - `window`
        - `.q_transform`
            - `frange`
            - `qrange`
        - `.gate`  
        To remove glitches via 'gating' method, which is padding basically center$\pm$width region
            - `tzero`  
            where to pad (center)
            - `tpad`  
            how many time to pad (width)

- `time`  
For converting usual `yymmdd` time to `gps time`

## `pycbc` (2.0.5) (`lalsuite` (7.11))
For generating waveforms, template generation, matched filtering, template bank generation - waveform template related stuffs

- `waveform`  
[td: time domain, fd: frequency domain]
    - `td_approximants`
    - `fd_approximants`  
    shows the available approximants
    - `get_td_waveform`
    - `get_fd_waveform`  
    Returns TimeSeries/FrequencySeries object of two polarization hplus and hcross;  
    parameters: `approximant`, `mass1`, `mass2`, `delta_t`/`delta_f`, `f_lower`, `distance` (see full documentation for other parameters)

- `psd`
    - `aLIGOZeroDetHighPower`
    - `welch`
    - `interpolate`
    - `inverse_spectrum_truncation`

- `noise`
    - `noise_from_psd`

- `types`
    - `TimeSeries`
        - `.crop`
        - `.psd`
        - `.cyclic_time_shift`
        - `.sample_times`
        - `.time_slice`
        - `.roll`   (shifts the whole data by some indices)
        - `.sample_frequencies`
        - `.to_frequencyseries`
        - `.to_timeseries`
        - `.start_time`
        - `.highpass_fir`
        - `.lowpass_fir`
        - `.time_slice`
        - `.whiten`
        - `.qtransform`
            - `logfsteps`
            - `qrange`
            - `frange`

- `catalog`
    - `Merger`  
    Gets merger data from [gwosc](https://gwosc.org).  
        - `.strain`  
        Returns TimeSeries object

- `frame`
    - `read_frame`  
    to read local file, expects `file`, `channel_name`.  
    Returns TimeSeries object

- `filter`
    - `resample_to_delta_t`
    - `highpass`
    - `matched_filter`  
    Returns TimeSeries object
    - `sigma`

- `vetoes`
    - `power_chisq`

- `events`
    - `ranking`
        - `newsnr`  
        Returns numpy array, not TimeSeries object

- `detector`
    - `Detector`
        - `light_travel_time_to_detctor`


## `bilby` (2.0.1)  
For Parameter estimation

- `core`
    - `likelihood`
        - `GaussianLikelihood`
    - `prior`
        - `Uniform`
        - `PowerLaw`
        - `PriorDict`
    - `result`
        - `Result`
            - `plot_corner`
            - `posterior`
            - `priors`
            - `sampler_kwargs`
            - `log_bayes_factor`
            - `log_evidence_err`
        - `read_in_result`  
        Reads from already generated and exported result, returns Result object
    
    - `utils`
        - `SamplesSummary`
            - `.median`
            - `.upper_relative_credible_interval`
            - `.lower_relative_credible_interval`
            - `.upper_absolute_credible_interval`
            - `.lower_absolute_credible_interval`

- `run_sampler`  
    Returns Result object

- `gw`
    - `conversion`
        - `convert_to_lal_binary_black_hole_parameters`
        - `generate_all_bbh_parameters`
    - `detector`
        - `get_empty_interferometer`  
        Returns Interferometer object
        - `interferometer`
            - `Interferometer`
                - `set_strain_data_from_gwpy_timeseries`
                - `strain_data`
                    - `roll_off`
                    - `frequency_mask`
                    - `frequency_array`
                    - `frequency_domain_strain`
                - `power_spectral_density`  
                This is a PowerSpectralDensity object
                - `maximum_frequency`
                - `minimum_frequency`
        - `psd`
            - `PowerSpectralDensity`
                - `frequency_array`
                - `asd_array`
        - `PowerSpectralDensity` (alias to psd.PowerSpectralDensity)
    - `WaveformGenerator`
    - `source`
        - `lal_binary_black_hole`
    - `likelihood`
        - `GravitationalWaveTransient`

## `dynesty` (2.1.1)


## `h5py`
To read `.hdf5` files

- `File` (alias to `_hl.files.File`)
- `._hl`
    - `files`
        - `File`
            - `keys`
                * `luminosity_distance_Mpc`: luminosity distance [Mpc]

                * `m1_detector_frame_Msun`: primary (larger) black hole mass (detector frame) [solar mass]

                * `m2_detector_frame_Msun`: secondary (smaller) black hole mass (detector frame) [solar mass]

                * `right_ascension`, `declination`: right ascension and declination of the source [rad].

                * `costheta_jn`: cosine of the angle between line of sight and total angular momentum vector of the system.

                * `spin1`, `costilt1`: primary (larger) black hole spin magnitude (dimensionless) and cosine of the zenith angle between the spin and the orbital angular momentum vector of the system.

                * `spin2`, `costilt2`: secondary (smaller) black hole spin magnitude (dimensionless) and cosine of the zenith angle between the spin and the orbital angular momentum vector of the system.

