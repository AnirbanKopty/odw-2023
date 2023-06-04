from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters, generate_all_bbh_parameters
from bilby.core.prior import Uniform, PowerLaw
import bilby
from gwpy.timeseries import TimeSeries
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

start_time = datetime.now()
# trigger times
trg_time, trg_nsnr, trg_mass = np.genfromtxt(
    "triggers_PE.dat", delimiter=',', comments='#', unpack=True)

# ? control parameter
indices = range(len(trg_time))

H1strain = TimeSeries.read("challenge3_2048hz.gwf", "H1:CHALLENGE3")
L1strain = TimeSeries.read("challenge3_2048hz.gwf", "L1:CHALLENGE3")

# analysis time window
post_trigger_duration = 2
duration = 4
analysis_start = trg_time + post_trigger_duration - duration
psd_duration = duration * 32
psd_start_time = analysis_start - psd_duration

# to write to a file the result
file = open('./results/result.txt', 'w')
string = "Estimated Parameter values"
file.write(string + '\n' + '(with 90% confidence interval)\n')
file.write('='*(len(string)+5) + '\n\n')

# Analysis part
for i in indices:

    # Strain
    idx1 = int(H1strain.sample_rate.value * analysis_start[i])
    idx2 = int(H1strain.sample_rate.value * (analysis_start[i] + duration))

    # Interferometers
    H1 = bilby.gw.detector.get_empty_interferometer("H1")
    L1 = bilby.gw.detector.get_empty_interferometer("L1")
    H1.set_strain_data_from_gwpy_timeseries(H1strain[idx1:idx2])
    L1.set_strain_data_from_gwpy_timeseries(L1strain[idx1:idx2])

    # PSD
    psd_idx1 = int(H1strain.sample_rate.value * psd_start_time[i])
    psd_idx2 = int(H1strain.sample_rate.value *
                   (psd_start_time[i] + psd_duration))

    H1_psd = H1strain[psd_idx1:psd_idx2]
    L1_psd = L1strain[psd_idx1:psd_idx2]

    psd_alpha = 2 * H1.strain_data.roll_off / duration
    H1_psd = H1_psd.psd(fftlength=duration, overlap=0,
                        window=("tukey", psd_alpha), method="median")
    L1_psd = L1_psd.psd(fftlength=duration, overlap=0,
                        window=("tukey", psd_alpha), method="median")

    H1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=H1_psd.frequencies.value, psd_array=H1_psd.value)
    L1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=L1_psd.frequencies.value, psd_array=L1_psd.value)

    H1.maximum_frequency = 1024
    L1.maximum_frequency = 1024

    # Parameter Estimation part
    prior = bilby.core.prior.PriorDict()
    prior['chirp_mass'] = Uniform(name='chirp_mass', minimum=8.5, maximum=43.5)
    prior['mass_ratio'] = Uniform(name='mass_ratio', minimum=0.5, maximum=1)
    prior['phase'] = Uniform(name="phase", minimum=0, maximum=2*np.pi)
    prior['geocent_time'] = Uniform(
        name="geocent_time", minimum=trg_time[i]-0.1, maximum=trg_time[i]+0.1)
    prior['a_1'] = 0.0
    prior['a_2'] = 0.0
    prior['tilt_1'] = 0.0
    prior['tilt_2'] = 0.0
    prior['phi_12'] = 0.0
    prior['phi_jl'] = 0.0
    # some random values, won't matter, we are after mass determination only
    prior['dec'] = -1.2232
    prior['ra'] = 2.19432
    prior['theta_jn'] = 1.89694
    prior['psi'] = 0.532268
    prior['luminosity_distance'] = PowerLaw(
        alpha=2, name='luminosity_distance', minimum=50, maximum=2000, unit='Mpc', latex_label='$d_L$')

    # First, put our "data" created above into a list of interferometers (the order is arbitrary)
    interferometers = [H1, L1]

    # Next create a dictionary of arguments which we pass into the LALSimulation waveform - we specify the waveform approximant here
    waveform_arguments = dict(
        waveform_approximant='IMRPhenomD', reference_frequency=100., catch_waveform_errors=True)

    # Next, create a waveform_generator object. This wraps up some of the jobs of converting between parameters etc
    waveform_generator = bilby.gw.WaveformGenerator(
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        waveform_arguments=waveform_arguments,
        parameter_conversion=convert_to_lal_binary_black_hole_parameters)

    # Finally, create our likelihood, passing in what is needed to get going
    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers, waveform_generator, priors=prior,
        time_marginalization=True, phase_marginalization=True, distance_marginalization=True)

    result_bayesian = bilby.run_sampler(
        likelihood, prior, sampler='dynesty', outdir=f'bayesian_output', 
        label=f"CHALLENGE4PE",
        conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,
        # <- Arguments are used to make things fast - not recommended for general use
        nlive=250, dlogz=1.,
        clean=True)
    # result_bayesian = bilby.core.result.read_in_result("bayesian_output/CHALLENGE4PE_result.json")


    # documenting results
    string = f'The Signal at {trg_time[i]} s with nsnr {trg_nsnr[i]} with mass {trg_mass[i]} Msun'
    file.write(string + '\n' + '-'*len(string) + '\n')

    result_bayesian.plot_corner(
        parameters=["mass_1", "mass_2"], prior=True, save=False)    # since, save=True closes the figure handle
    plt.savefig(f'./results/corner{i}.png')

    mass = bilby.core.utils.SamplesSummary(
        samples=result_bayesian.posterior['mass_1'], average='median')
    file.write(
        f'mass_1 = {mass.median:.4f} +{mass.upper_relative_credible_interval:.4f} {mass.lower_relative_credible_interval:.4f} Msun\n')

    mass = bilby.core.utils.SamplesSummary(
        samples=result_bayesian.posterior['mass_2'], average='median')
    file.write(
        f'mass_2 = {mass.median:.4f} +{mass.upper_relative_credible_interval:.4f} {mass.lower_relative_credible_interval:.4f} Msun\n')
    file.write('\n')


end_time = datetime.now()

file.write(f"\n\n\n\n\nThe time taken for the script to be executed {end_time - start_time}")
file.close()