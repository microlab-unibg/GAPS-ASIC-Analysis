import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

import os

from functions.utils import adjacent_values

""" === PEDESTAL ===
Dispersion of output value of ADC with $Q_{inj} = 0C$.

### Measures
- Histogram (ADC values, Count)
- Mean, std dev by $\tau_{peak}$
- Mean, std dev by $\#CHANNEL$
- Mean, std dev for the chip
- Normality check (not normal)"""


def pedestal(isProcessed, enablePlot, singleMeasures, peaking_time, input_path, out_path_current, n_ch, n_tau, n_fthr):
    list_tau = [peaking_time] if singleMeasures else np.arange(n_tau)

    for j in list_tau:
        fname = input_path + "Pedestals_tau{}.dat".format(j)
        break

    try:
        val = np.loadtxt(fname, comments='#', usecols=(4), unpack=True)
    except OSError:
        print('No Pedestal file was found. Please check that pedestal file name meets "Pedestals_tauX.dat" pattern.')
        exit(1)

    pedestal_events = int(len(val) / n_ch)
    m_values = np.zeros((n_tau * pedestal_events, n_ch))

    # Processed files
    sep = os.path.sep
    processed_path = out_path_current + 'Processed' + sep
    if (not(os.path.exists(processed_path))):
        os.makedirs(processed_path)
    file_mu = processed_path + 'mu.dat'
    file_sigma = processed_path + 'sigma.dat'
    file_outliers = processed_path + 'outliers.dat'

    # Processed
    m_mu = np.zeros((n_tau, n_ch), dtype=float)
    m_sigma = np.zeros((n_tau, n_ch), dtype=float)
    m_outliers = np.zeros((n_tau, n_ch), dtype=float)

    if (isProcessed):
        try:
            m_mu = np.loadtxt(file_mu, dtype='float', delimiter='\t')
            m_sigma = np.loadtxt(file_sigma, dtype='float', delimiter='\t')
            m_outliers = np.loadtxt(file_outliers, dtype='float', delimiter='\t')
        except OSError:
            print('There are no processed files or some are missing for this analysis')
            exit(1)

        for j in list_tau:
            base_idx = j * pedestal_events
            # Getting data
            try:
                fname = input_path + "Pedestals_tau{}.dat".format(j)
                val = np.loadtxt(fname, dtype='int', comments='#', usecols=(4), unpack=True)
                m_values[base_idx:base_idx + pedestal_events, :] = val.reshape((pedestal_events, n_ch))
            except OSError:
                print('File Pedestals_tau{}.dat not found! Ignoring this analysis...'.format(j))
                pass
    else:
        """=== PROCESSING ==="""
        print('Processing input files...')
        for j in list_tau:
            base_idx = j * pedestal_events
            # Getting data
            try:
                fname = input_path + "Pedestals_tau{}.dat".format(j)
                typea, val = np.loadtxt(fname, comments='#', usecols=(2, 4), unpack=True)
                print('File Pedestal_tau{}.dat found!'.format(j))
                val[np.where((typea == 1) | (typea == 11))] = np.nan
                m_values[base_idx:base_idx + pedestal_events, :] = val.reshape((pedestal_events, n_ch))
            except OSError:
                print('File Pedestal_tau{}.dat not found! Ignoring this analysis...'.format(j))
                pass

        with open(out_path_current + 'Pedestal.dat', 'w') as f:
            for i in range(n_ch):
                print('CH {}:'.format(i))
                f.write('CH {}:\n'.format(i))
                values_bp = []
                num_nans = 0
                for j in list_tau:
                    # Mean, std dev by tau
                    base_idx = j * pedestal_events
                    this_vals = np.copy(m_values[base_idx:base_idx + pedestal_events, i])
                    num_nans += len(this_vals[np.isnan(this_vals)])

                    (mu, sigma) = np.nanmean(this_vals), np.nanstd(this_vals)
                    m_mu[j, i] = mu
                    m_sigma[j, i] = sigma
                    f.write('Tau {}: mean={}, sigma={:.2f}\n'.format(j, mu, sigma))
                    print('Tau', j, ': mean', mu, ', sigma', '{:.2f}'.format(sigma))

                    # Normality check (ricontrolla)
                    """values_norm = (values[j, idx] - np.full(np.size(values[j, idx]), mu))/(sigma)
                    k2, p = stats.normaltest(values_norm)
                    sm.qqplot(values_norm, line = '45')
                    print("Normal-like" if p > alpha else "Not normal")
                    """

                    values_bp.append(this_vals[~np.isnan(this_vals)].flatten())

                # Outliers
                fliers = plt.boxplot(values_bp)['fliers']
                for j in range(len(list_tau)):
                    # the y and x positions of the fliers
                    yfliers = fliers[j].get_ydata()
                    m_outliers[j, i] = len(yfliers)

                plt.close()

                # Mean, std dev by channel
                mu_ch, sigma_ch = np.mean(m_mu, axis=0), np.mean(m_sigma, axis=0)
                f.write('Channel {}: mean={:.2f}, sigma={:.2f}\n'.format(i, mu_ch[i], sigma_ch[i]))
                print('Channel', i, ': mean', '{:.2f}'.format(mu_ch[i]), ', sigma', '{:.2f}'.format(sigma_ch[i]))

                f.write('Number of error or empty events: {}\n'.format(num_nans))
                print('Number of error or empty events: {}'.format(num_nans))

            # Mean, std dev for the chip
            mu_chip, sigma_chip = np.mean(mu_ch), np.mean(sigma_ch)
            f.write('Chip: mean={:.2f}, sigma={:.2f}\n'.format(mu_chip, sigma_chip))
            print('Chip: mean', '{:.2f}'.format(mu_chip), ', sigma', '{:.2f}'.format(sigma_chip))

        # Save processed data
        header = 'Peak times (rows) and Channels (cols)\n'
        header = header + '\t'.join(['ch{}'.format(i) for i in range(n_ch)])
        np.savetxt(file_mu, m_mu, delimiter='\t', header=header)
        np.savetxt(file_sigma, m_sigma, delimiter='\t', header=header)
        np.savetxt(file_outliers, m_outliers, delimiter='\t', header=header)

    if (enablePlot):
        """=== PLOTTING ==="""
        sep = os.path.sep
        out_path_old = out_path_current + 'Pedestal_Histogram' + sep
        if (not(os.path.exists(out_path_old))):
            os.makedirs(out_path_old)
        out_path_violin = out_path_current + 'Pedestal_Violin' + sep
        if (not(os.path.exists(out_path_violin))):
            os.makedirs(out_path_violin)
        out_path_boxplot = out_path_current + 'Pedestal_Boxplot' + sep
        if (not(os.path.exists(out_path_boxplot))):
            os.makedirs(out_path_boxplot)
        print('Plotting...')
        for i in range(n_ch):
            print('CH {}:'.format(i))
            out_path_single = out_path_old + 'Channel_' + str(i) + sep
            if (not(os.path.exists(out_path_single))):
                os.makedirs(out_path_single)

            values_bp = []
            bp_labels = []
            fig, ax1 = plt.subplots(figsize=(10, 8))
            for j in list_tau:
                base_idx = j * pedestal_events
                print('  Tau {}:'.format(j))
                fig_s, ax_s = plt.subplots(figsize=(10, 8))

                # Histogram
                this_vals = m_values[base_idx:base_idx + pedestal_events, i]
                this_vals = this_vals[~np.isnan(this_vals)]
                values_bp.append(this_vals)
                print('    Histogram')
                xmin = int(np.nanmin(this_vals)) - 10
                xmax = int(np.nanmax(this_vals)) + 10

                bp_labels.append("$\\tau_{}$\n$\mu={:.2f}$\n$\sigma={:.2f}$".format(j, m_mu[j, i], m_sigma[j, i]))

                h, b, patches = ax1.hist(this_vals, range=(xmin, xmax), bins=(xmax - xmin + 1), histtype="bar",
                                         alpha=0.5, edgecolor='black', linewidth=1.2,
                                         label="$\\tau_{}\quad(\mu={:.2f},\sigma={:.2f}$)".format(j, m_mu[j, i], m_sigma[j, i]))

                ax1.set_ylabel('Counts')

                h_s, b_s, patches_s = ax_s.hist(this_vals, range=(xmin, xmax), bins=(xmax - xmin + 1), histtype="bar",
                                                fc=(0, 0, 1, 0.5), edgecolor='black', linewidth=1.2,
                                                label="$\\tau_{}\quad(\mu={:.2f},\sigma={:.2f}$)".format(j, m_mu[j, i], m_sigma[j, i]))
                ax_s.set_ylabel('Counts')
                ax_s.set_xlabel("Channel_out [ADC code]")
                chartBox = ax_s.get_position()
                ax_s.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.8, chartBox.height])
                ax_s.legend(loc=7, bbox_to_anchor=(1.385, 0.5), borderaxespad=0, frameon=True, ncol=1)
                plt.grid(True)
                plt.title("Pedestal of channel #{}".format(i))
                plt.savefig(out_path_single + 'Pedestal_ch' + str(i) + '_tau' + str(j) + '.svg', format='svg', bbox_inches="tight")
                plt.close()

            # Outliers & Box Plot
            print('  Boxplot')
            fig_bp, ax1_bp = plt.subplots(figsize=(11, 8))
            ax1_bp.set_title('Boxplot for channel #{}'.format(i))
            fliers = ax1_bp.boxplot(values_bp, labels=bp_labels)['fliers']
            for k in range(len(fliers)):
                # the y and x positions of the fliers
                yfliers = fliers[k].get_ydata()
                bp_labels[k] = bp_labels[k] + "\n#Outliers={}".format(len(yfliers))

            ax1_bp.boxplot(values_bp, labels=bp_labels)
            ax1_bp.set_xlabel('$\\tau_p$')
            ax1_bp.set_ylabel('Channel_out [ADC code]')

            plt.savefig(out_path_boxplot + 'Pedestal_ch' + str(i) + '_boxplot.svg', format='svg', bbox_inches="tight")
            plt.close()

            # Violin plot
            # courtesy of https://matplotlib.org/gallery/statistics/customized_violin.html#sphx-glr-gallery-statistics-customized-violin-py)
            print('  Violinplot')
            fig_vp, ax1_vp = plt.subplots(figsize=(11, 8))
            ax1_vp.set_title('Violin plot for channel #{}'.format(i))

            ax1_vp.violinplot(values_bp)

            quartile1, medians, quartile3 = np.percentile(values_bp, [25, 50, 75], axis=1)
            whiskers = np.array([
                adjacent_values(sorted_array, q1, q3)
                for sorted_array, q1, q3 in zip(values_bp, quartile1, quartile3)])
            whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]

            inds = np.arange(1, len(medians) + 1)
            ax1_vp.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
            ax1_vp.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
            ax1_vp.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)

            ax1_vp.set_xlabel('$\\tau_p$')
            ax1_vp.set_xticks(np.arange(1, len(bp_labels) + 1))
            ax1_vp.set_xticklabels(bp_labels)
            ax1_vp.set_ylabel('Channel_out [ADC code]')

            plt.savefig(out_path_violin + 'Pedestal_ch' + str(i) + '_violinplot.svg', format='svg', bbox_inches="tight")
            plt.close()

            # Showing the plot
            print('  Channel Histogram')
            mu_ch, sigma_ch = np.mean(m_mu[:, i]), np.mean(m_sigma[:, i])
            x = np.linspace(mu_ch - 3 * sigma_ch, mu_ch + 3 * sigma_ch, 100)
            ax2 = ax1.twinx()
            ax2.plot(x, stats.norm.pdf(x, mu_ch, sigma_ch), color='b', linewidth=1, label=r"N$(\mu={:.2f},\sigma={:.2f})$".format(mu_ch, sigma_ch))
            ax2.set_ylabel("Norm")
            ax1.set_xlabel("Channel_out [ADC code]")
            chartBox = ax1.get_position()
            ax1.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.85, chartBox.height])
            ax1.legend(loc=7, bbox_to_anchor=(1.55, 0.6), borderaxespad=0, frameon=True, ncol=1)
            # ax1.legend()
            chartBox = ax2.get_position()
            ax2.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.85, chartBox.height])
            ax2.legend(loc=7, bbox_to_anchor=(1.518, 0.425), borderaxespad=0, frameon=True, ncol=1)
            #ax2.legend(loc='upper left')
            plt.title("Pedestal of channel #{}".format(i))
            plt.grid(True)
            # plt.show()
            if not(singleMeasures):
                plt.savefig(out_path_single + 'Pedestal_ch' + str(i) + '_histogram.svg', format='svg', bbox_inches="tight")
            plt.close()

    return m_sigma
