import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as plt_colours

""" === WAVEFORM ===
Shaper output waveform

### Measures
- Waveform
- $\tau_{peak}$ corresponding to max ADC code
- Mean, std dev by $\tau_{peak}$ using first `n_sample` samples
- Mean, std dev by $\tau_{peak}$ using mean and std dev of max ADC code
- Mean, std dev for the chip using both methods"""


def waveform(enablePlot, input_path, out_path_current, n_ch, n_tau, n_fthr):
    # Number of highest samples
    n_sample = 32

    colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

    # Getting number of samples
    fname = input_path + "WaveformScan_fast_tau0.dat"
    try:
        dta = np.loadtxt(fname, dtype='int', comments='#', usecols=(0), unpack=True)
    except OSError:
        pass

    waveform_sample = len(np.unique(dta))
    m_mu = np.zeros((n_tau * waveform_sample, n_ch))
    m_sigma = np.zeros((n_tau * waveform_sample, n_ch))
    m_t_peak_mu = np.zeros((n_tau, n_ch))
    m_maximum_out_mu = np.zeros((n_tau, n_ch))
    m_maximum_out_sigma = np.zeros((n_tau, n_ch))
    m_t_peak_mu_first = np.zeros((n_tau, n_ch))
    m_t_peak_sigma_first = np.zeros((n_tau, n_ch))
    m_t_peak_mu_second = np.zeros((n_tau, n_ch))
    m_t_peak_sigma_second = np.zeros((n_tau, n_ch))

    # Output file format
    header = 'Peak times (rows) and Channels (cols)\n'
    header = header + 'tp\\ch\t' + '\t'.join(['ch{}'.format(i) for i in range(n_ch)])
    data = np.zeros((n_tau, 1 + n_ch + 2))
    data[0:n_tau, 0] = [str(j) for j in range(n_tau)]
    data_max = np.zeros((n_tau, 1 + n_ch + 2))
    data_max[0:n_tau, 0] = [str(j) for j in range(n_tau)]
    data1 = np.zeros((n_tau, 1 + n_ch + 2))
    data1[0:n_tau, 0] = [str(j) for j in range(n_tau)]
    data2 = np.zeros((n_tau, 1 + n_ch + 2))
    data2[0:n_tau, 0] = [str(j) for j in range(n_tau)]
    # Output path
    file_waveform = out_path_current + 'waveform.dat'
    file_waveform_max = out_path_current + 'waveform_max.dat'
    file_waveform_first = out_path_current + 'waveform_first.dat'
    file_waveform_second = out_path_current + 'waveform_second.dat'

    """=== PROCESSING ==="""
    for i in range(n_ch):
        plt.figure(figsize=(8, 6))
        ax = plt.subplot(111)
        ax.set_ylim(0, 1000)
        print('CH {}:'.format(i))
        for j in range(n_tau):
            print('  Tau {}:'.format(j))
            base_idx = j * waveform_sample

            # Getting data
            try:
                fname = input_path + "WaveformScan_fast_tau" + str(j) + ".dat"
                dta, daca, typea, ch, val = np.loadtxt(fname, dtype='int', comments='#', usecols=(0, 1, 2, 3, 4), unpack=True)
                print('File WaveformScan_fast_tau{}.dat found!'.format(j))
                dt = np.unique(dta)
                dac = np.unique(daca)
            except OSError:
                print('File WaveformScan_fast_tau{}.dat not found! Ignoring...'.format(j))
                pass

            idx = (ch == i).nonzero()[0]
            val[np.where((typea == 1) | (typea == 11))] = np.nan

            # Mean, std dev
            for k in range(len(dt)):
                jdx = (dta[idx] == dt[k]).nonzero()[0]
                mu, sigma = np.mean(val[idx[jdx]]), np.std(val[idx[jdx]])
                m_mu[base_idx + k, i] = mu
                m_sigma[base_idx + k, i] = sigma

            # Mean, std dev by tau w/ max samples (48 is a conversion factor)
            this_mu = np.copy(m_mu[base_idx:base_idx + waveform_sample, i])
            this_sigma = np.copy(m_sigma[base_idx:base_idx + waveform_sample, i])
            maximum_idx = np.argmax(this_mu)
            maximum_tpeak = (maximum_idx + 1) / 48
            maximum_out_mu = this_mu[maximum_idx]
            maximum_out_sigma = this_sigma[maximum_idx]
            m_t_peak_mu[j, i] = maximum_tpeak
            m_maximum_out_mu[j, i] = maximum_out_mu
            m_maximum_out_sigma[j, i] = maximum_out_sigma
            print('    Max t_peak: {:.2f}'.format(maximum_tpeak))

            # order by peak
            dta = dta / 48
            val_sorted = sorted(val, reverse=True)
            idxs = np.where(val >= val_sorted[n_sample - 1])
            t_peak_mu, t_peak_sigma = np.mean(dta[idxs]), np.std(dta[idxs])
            m_t_peak_mu_first[j, i] = t_peak_mu
            m_t_peak_sigma_first[j, i] = t_peak_sigma
            print('      First: t_peak mu: {:.2f}'.format(t_peak_mu), 't_peak sigma: {:.2f}'.format(t_peak_sigma))

            # Mean, std dev by tau w/ ADC code mean, std dev
            idxs = [l for l in range(len(val))
                    if (val[l] >= maximum_out_mu - 3 * maximum_out_sigma) & (val[l] <= maximum_out_mu + 3 * maximum_out_sigma)]
            t_peak_mu, t_peak_sigma = np.mean(dta[idxs]), np.std(dta[idxs])
            m_t_peak_mu_second[j, i] = t_peak_mu
            m_t_peak_sigma_second[j, i] = t_peak_sigma
            print('      Second: t_peak mu: {:.2f}'.format(t_peak_mu), 't_peak sigma: {:.2f}'.format(t_peak_sigma))

    # Saving data for file writing
    data[0:n_tau, 1:n_ch + 1] = m_t_peak_mu
    data_max[0:n_tau, 1:n_ch + 1] = m_maximum_out_mu
    data1[0:n_tau, 1:n_ch + 1] = m_t_peak_mu_first
    data2[0:n_tau, 1:n_ch + 1] = m_t_peak_mu_second

    # Mean, std dev for the chip for both methods
    print('Chip')
    for j in range(n_tau):
        print('  tau', j)
        t_peak_mu_chip = np.mean(m_t_peak_mu[j, :])
        print('    raw: t_peak mu:', '{:.2f}'.format(t_peak_mu_chip))
        maximum_mu_chip, maximum_sigma_chip = np.mean(m_maximum_out_mu[j, :]), np.average(m_maximum_out_sigma[j, :], weights=m_maximum_out_mu[j, :])
        print('    maximum: mu:', '{:.2f}'.format(maximum_mu_chip), 'sigma:', '{:.2f}'.format(maximum_sigma_chip))
        t_peak_mu_first_chip, t_peak_sigma_first_chip = np.mean(m_t_peak_mu_first[j, :]), np.average(m_t_peak_sigma_first[j, :], weights=m_t_peak_mu_first[j, :])
        print('    First: t_peak mu:', '{:.2f}'.format(t_peak_mu_first_chip), 't_peak sigma:', '{:.2f}'.format(t_peak_sigma_first_chip))
        t_peak_mu_second_chip, t_peak_sigma_second_chip = np.mean(m_t_peak_mu_second[j, :]), np.average(m_t_peak_sigma_second[j, :], weights=m_t_peak_mu_second[j, :])
        print('    Second: t_peak mu:', '{:.2f}'.format(t_peak_mu_second_chip), 't_peak sigma:', '{:.2f}'.format(t_peak_sigma_second_chip))

        data[j, n_ch + 1] = t_peak_mu_chip
        data_max[j, n_ch + 1], data_max[j, n_ch + 2] = maximum_mu_chip, maximum_sigma_chip
        data1[j, n_ch + 1], data1[j, n_ch + 2] = t_peak_mu_first_chip, t_peak_sigma_first_chip
        data2[j, n_ch + 1], data2[j, n_ch + 2] = t_peak_mu_second_chip, t_peak_sigma_second_chip

    # File writing
    np.savetxt(file_waveform, data, delimiter='\t', header=header)
    np.savetxt(file_waveform_max, data_max, delimiter='\t', header=header)
    np.savetxt(file_waveform_first, data1, delimiter='\t', header=header)
    np.savetxt(file_waveform_second, data2, delimiter='\t', header=header)

    """=== PLOTTING ==="""
    if (enablePlot):
        print('Plotting tau for channel...')
        idx_max_t_peak = np.argmax(m_mu)
        upper_limit = (m_mu.flatten())[idx_max_t_peak] + np.multiply(3, (m_sigma.flatten())[idx_max_t_peak])
        # Plot tau for each channel
        for i in range(n_ch):
            plt.figure(figsize=(8, 6))
            ax = plt.subplot(111)
            ax.set_ylim(0, upper_limit)
            print('  Channel {}:'.format(i))
            for j in range(n_tau):
                base_idx = j * waveform_sample

                ax.plot((dt + 1) / 48, m_mu[base_idx:base_idx + waveform_sample, i], label="$\\tau_{}, {}={:.2f}\mu s$".format(j, '\\tau_{p}', maximum_tpeak), c=colours[j])

                # Confidence intervals
                y_conf_up = m_mu[base_idx:base_idx + waveform_sample, i] + np.multiply(3, m_sigma[base_idx:base_idx + waveform_sample, i])
                y_conf_down = m_mu[base_idx:base_idx + waveform_sample, i] - np.multiply(3, m_sigma[base_idx:base_idx + waveform_sample, i])

                colour = plt_colours.to_rgba(colours[j])
                colour = list(colour)
                colour[3] = colour[3] - 0.7
                ax.plot((dt + 1) / 48, y_conf_up, c=colour)
                ax.plot((dt + 1) / 48, y_conf_down, c=colour)
                ax.fill_between((dt + 1) / 48, y_conf_down, y_conf_up, color=colour)

            # Show the plot
            plt.xlabel("t [$\\mu$s]")
            plt.ylabel("Channel_out [ADC code]")
            plt.title("Waveform of channel #{} ($DAC_{{inj}}$ code = {})".format(i, dac[0]))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.85, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.3, 0.5), borderaxespad=0, frameon=True, ncol=1)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'Waveform_ch' + str(i) + '.svg', format='svg', bbox_inches="tight")
            plt.close()


        # Plot channels for each tau
        print('Plotting channels for tau...')
        for j in range(n_tau):
            print('  Tau {}'.format(j))
            base_idx = j * waveform_sample
            plt.figure(figsize=(12, 11))
            ax = plt.subplot(111)
            ax.set_ylim(0, upper_limit)
            for i in range(n_ch):
                ax.plot((dt + 1) / 48, m_mu[base_idx:base_idx + waveform_sample, i], label="CH #{}, $\\tau_{{p}}={:.2f}\mu s$".format(i, data[j, i + 1]))

            plt.xlabel("t [$\\mu$s]")
            plt.ylabel("Channel_out [ADC code]")
            plt.title("Waveform for $\\tau_{}$ ($DAC_{{inj}}$ code = {})".format(j, dac[0]))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.9, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.24, 0.5), borderaxespad=0, frameon=True, ncol=1)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'Waveform_tau' + str(j) + '_allch.svg', format='svg', bbox_inches="tight")
            plt.close()
