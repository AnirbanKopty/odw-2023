# Detectors
- `'G1'` - GEO600
- `'H1'` - LIGO-Hanford
- `'L1'` - LIGO-Livingston
- `'V1'` - (Advanced) Virgo
- `'K1'` - KAGRA


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
        - `fetch_open_data`  
        TimeSeries object that has 1D data of GW signal for a particular time segment for a particular detector.  
        Creates a TimeSeries object

        Has several methods
        - `.plot`
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

- `time`  
For converting usual `yymmdd` time to `gps time`

## `pycbc` (2.0.5) (`lalsuite` (7.11))
For generating waveforms, template generation, matched filtering, template bank generation - waveform template related stuffs

- `waveform`
    - `td_approximants`
    - `fd_approximants`
    - `get_td_waveform`
    - `get_fd_waveform`  
    Returns TimeSeries object

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
        - `.strain`  
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

- `detector`
    - `Detector`
        - `light_travel_time_to_detctor`