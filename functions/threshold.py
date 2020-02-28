import numpy as np
from scipy import optimize
import math

import matplotlib.pyplot as plt
from matplotlib import ticker

import os

from functions.utils import f_VTP, write_default

"""=== THRESHOLD SCAN ===
Varying DAC voltage, plot how many times we trigger an event

### Measures
- Threshold with fixed tau and variable fthr (not implemented, at the end of the file)
- Threshold with fthr fixed and variable peaking time
- Threshold position (for each fine thr level and for each channel)
- Dispersion (for each channel)"""


def threshold(enablePlot, input_path, out_path_current, n_ch, n_tau, n_fthr):
    # with open(out_path_current + 'ThresholdScan.dat','w') as f:
    for i in range(n_ch):
        # f.write('CH #{}'.format(i))
        print('CH #{}'.format(i))
        for fthr in range(n_fthr):
            plt.figure(figsize=(9, 6))
            ax = plt.subplot(111)
            #f.write('  fthr: {}'.format(fthr))
            print('  fthr: {}'.format(fthr))
            vec_th = []
            some_nans = 0

            for j in range(n_tau):
                ok = 0
                # Getting data
                try:
                    fname = input_path + "ThresholdScan_fthr" + str(fthr) + "_tau" + str(j) + ".dat"
                    threshold, daca, events, triggered, ch = np.loadtxt(fname, dtype='int', comments='#', usecols=(0, 1, 2, 3, 4), unpack=True)
                except OSError:
                    pass
                thr = np.unique(threshold)
                ev = np.unique(events)
                dac = np.unique(daca)
                y = np.zeros(len(thr))
                ok = 1

                # Plotting data
                for k in range(len(thr)):
                    idx = (ch == i).nonzero()[0]
                    jdx = (threshold[idx] == thr[k]).nonzero()[0]
                    y[k] = 100 * triggered[idx[jdx]] / ev

                if (enablePlot):
                    ax.plot(thr, y, label="$\\tau_{}$".format(j))

            if (enablePlot):
                plt.xlabel("Discriminator threshold [$DAC_{thr}$ code]")
                plt.ylabel("Efficiency [%]")
                plt.title("Threshold scan of channel #{} with $DAC_{{inj}}$ code = {}".format(i, dac[0]))
                chartBox = ax.get_position()
                ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.95, chartBox.height])
                ax.legend(loc=7, bbox_to_anchor=(1.135, 0.5), borderaxespad=0, frameon=True, ncol=1)
                #plt.legend(loc = "lower right")
                plt.grid(True)
                # plt.show()
                plt.savefig(out_path_current + 'ThresholdScan_fthr{}_ch{}.svg'.format(fthr, i), format='svg', bbox_inches="tight")
                plt.close()


"""=== FINDING THE MINIMUM AND MAXIMUM OF FITTED VALUE ==="""


def find_min_max_threshold():
    print('Finding maximum and minimum fitting value for graph plotting...')

    plt_fit_xmin = 200
    plt_fit_xmax = 250

    return plt_fit_xmax, plt_fit_xmin


"""=== FITTING THRESHOLD SCAN FUNCTION ===
The error function is used to determinate the threshold voltage of the comparator and the std dev of the noise in the DAC
"""


def threshold_fitting_tau(isProcessed, enablePlot, input_path, out_path_current, n_ch, n_tau, n_fthr, plt_fit_xmax, plt_fit_xmin):
    # Processed files
    sep = os.path.sep
    processed_path = out_path_current + 'Processed' + sep
    if (not(os.path.exists(processed_path))):
        os.makedirs(processed_path)
    file_a_TAU = processed_path + 'a_tau.dat'
    file_b_TAU = processed_path + 'b_tau.dat'

    a = np.zeros((n_tau, n_ch * n_fthr), dtype=float)
    b = np.zeros((n_tau, n_ch * n_fthr), dtype=float)
    dac = []

    if (isProcessed):
        a = np.loadtxt(file_a_TAU, delimiter='\t', dtype=float)
        b = np.loadtxt(file_b_TAU, delimiter='\t', dtype=float)
        try:
            fname = input_path + "ThresholdScan_fthr0_tau0.dat"
            dac = np.loadtxt(fname, dtype='int', comments='#', usecols=(1), unpack=True)
        except OSError:
            print('There are no processed files or some are missing for this analysis')
            exit(1)
    else:
        """=== PROCESSING ==="""
        print('Processing input files...')
        # with open(out_path_current + 'ThresholdScan.dat','w') as f:
        for i in range(n_ch):
            # f.write('CH #{}\n'.format(i))
            print('CH {}:'.format(i))
            for fthr in range(n_fthr):
                #f.write('  fthr: {}\n'.format(fthr))
                print('  fthr: {}'.format(fthr))
                for j in range(n_tau):
                    #f.write('    tau: {}\n'.format(j))
                    print('    tau: {}'.format(j))
                    ok = 0
                    # Getting data
                    try:
                        fname = input_path + "ThresholdScan_fthr" + str(fthr) + "_tau" + str(j) + ".dat"
                        threshold, dac, events, triggered, ch = np.loadtxt(fname, dtype='int', comments='#', usecols=(0, 1, 2, 3, 4), unpack=True)
                    except OSError:
                        pass
                    thr = np.unique(threshold)
                    ev = np.unique(events)
                    y = np.zeros(len(thr))
                    ok = 1

                    # Plotting data
                    for k in range(len(thr)):
                        idx = (ch == i).nonzero()[0]
                        jdx = (threshold[idx] == thr[k]).nonzero()[0]
                        y[k] = triggered[idx[jdx]] / ev

                    # fitting function
                    b_init = np.random.normal(20, 1.5)
                    a_init = np.random.normal(250, 5)
                    x = thr
                    param = np.array([a_init, b_init])

                    try:
                        popt, pcov = optimize.curve_fit(f_VTP, x, y, param)  # bounds = (lower_bounds, upper_bounds))
                    except RuntimeError:
                        #f.write('    Not fitted for tau{}\n'.format(j))
                        print('    Not fitted for tau{}'.format(j))

                    # print(popt)
                    a[j, i * n_fthr + fthr] = popt[0]
                    b[j, i * n_fthr + fthr] = popt[1]

                    # Finding the threshold
                    #f.write('    Threshold for tau{}={:.2f}, b={:.2f}\n'.format(j, a_fthr[j],b_fthr[j]))
                    print('    Threshold for tau{}={:.2f}, b={:.2f}'.format(j, popt[0], popt[1]))

            #f.write('  Threshold for channel {}={:.2f}, sigma={:.2f}, dispersion={:.2f}\n'.format(i, a_mu, a_sigma, 3*2*a_sigma))
            a_mu, a_sigma = np.mean(a[:, i * n_fthr:i * n_fthr + n_fthr]), np.std(a[:, i * n_fthr:i * n_fthr + n_fthr])
            print('  Threshold for CH {}:'.format(i),
                  'mu={:.2f}'.format(a_mu),
                  ', sigma={:.2f}'.format(a_sigma),
                  ', dispersion:{:.2f}'.format(3 * 2 * a_sigma))

        np.savetxt(file_a_TAU, a, delimiter='\t')
        np.savetxt(file_b_TAU, b, delimiter='\t')

    if (enablePlot):
        """=== PLOTTING ==="""
        print('Plotting...')
        for i in range(n_ch):
            print('CH {}:'.format(i))
            for fthr in range(n_fthr):
                plt.figure(figsize=(9, 6))
                ax = plt.subplot(111)
                #plt.xlim(np.floor(a_min) -5, np.floor(a_max) + 5)
                print('  fthr: {}'.format(fthr))
                for j in range(n_tau):
                    x = np.linspace(100, 300, 500)
                    ax.plot(x, 100 * f_VTP(x, a[j, i * n_fthr + fthr], b[j, i * n_fthr + fthr]), label="$\\tau_{}$, a={:.2f}, b={:.2f}".format(j, a[j, i * n_fthr + fthr], b[j, i * n_fthr + fthr]))

                a_mu, a_sigma = np.mean(a[:, i * n_fthr:i * n_fthr + n_fthr]), np.std(a[:, i * n_fthr:i * n_fthr + n_fthr])
                ax.errorbar(a_mu, 50, xerr=np.multiply(3, a_sigma), xlolims=True, xuplims=True, c='b')
                ax.plot(a_mu, 50, c='b', marker='o', linewidth=2)

                plt.xlim(plt_fit_xmin, plt_fit_xmax)
                plt.xlabel("Discriminator threshold [$DAC_{thr}$ code]")
                plt.ylabel("Efficiency [%]")
                bin_fthr = format(fthr, "03b")
                plt.title("Threshold of channel #{} with $DAC_{{inj}}$ = {} and $fthr = {}$".format(i, dac[0], bin_fthr))
                chartBox = ax.get_position()
                ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.8, chartBox.height])
                ax.legend(loc=7, bbox_to_anchor=(1.39, 0.5), borderaxespad=0, frameon=True, ncol=1)
                #plt.legend(loc = "lower right")
                plt.grid(True)
                # plt.show()
                plt.savefig(out_path_current + 'ThresholdScan_fthr{}_ch{}.svg'.format(fthr, i), format='svg', bbox_inches="tight")
                plt.close()

        # Plot for a (fitted voltage threshold) and b (standard deviation noise)"""
        # Plot for b
        for fthr in range(n_fthr):

            # Plot for b (x -> channels)
            plt.figure(figsize=(8, 6))
            plt.ylim(np.amin(b) - 2, np.amax(b) + 2)
            plt.xlim(-1, 32)
            ax = plt.gca()
            yticks = ax.xaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)

            x = np.linspace(0, 31, 32, endpoint=True)
            marker = ['-*', '-x', '-+', '-D', '-o', '-^', '-p', '-v']
            for j in range(n_tau):
                ax.plot(x, b[j, fthr::n_fthr], marker[j], markersize=5, linewidth=0.7, label="$\\tau_{}$".format(j))

            bin_fthr = format(fthr, "03b")
            plt.xlabel("Channel #")
            plt.ylabel("rms [$DAC_{thr}$ code]")
            plt.title("Noise standard deviation with fthr = {}".format(bin_fthr))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.9, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.19, 0.5), borderaxespad=0, frameon=True, ncol=1)
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.tick_params(axis='x', which='minor', direction='out', bottom=True, length=3)

            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'B_fthr{}_CH.svg'.format(fthr), format='svg', bbox_inches="tight")
            plt.close()

            # Plot for a
            plt.figure(figsize=(8, 8))
            ax = plt.subplot(111)
            plt.ylim(np.amin(a) - 5, np.amax(a) + 5)
            plt.xlim(-1, 32)
            plt.grid(True)

            for j in range(n_tau):
                ax.plot(x, a[j, fthr::n_fthr], marker[j], markersize=5, linewidth=0.7, label="$\\tau_{}$".format(j))

            plt.xlabel("Channel #")
            plt.ylabel("Discriminator threshold [$DAC_{thr}$ code]")
            plt.title("Threshold fthr = {}".format(bin_fthr))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.9, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.19, 0.5), borderaxespad=0, frameon=True, ncol=1)
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.tick_params(axis='x', which='minor', direction='out', bottom=True, length=3)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'A_fthr{}.svg'.format(fthr), format='svg', bbox_inches="tight")
            plt.close()

            # Plot for b (x -> tau)
            plt.figure(figsize=(10, 10))
            plt.ylim(np.amin(b) - 2, np.amax(b) + 2)
            plt.xlim(0, 8)
            ax = plt.gca()
            yticks = ax.xaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)

            x = np.linspace(0, 7, 8, endpoint=True)
            #marker = ['-*', '-x', '-+', '-D', '-o', '-^', '-p', '-v']
            #print("fthr: {}".format(fthr))
            for i in range(n_ch):
                    #print("Channel: {}".format(j))
                    #print( b[(n_tau*n_fthr*j)+fthr*n_tau:(n_tau*n_fthr*j)+fthr*n_tau+n_tau])
                ax.plot(x, b[:, i * n_fthr + fthr], label="CH #{}".format(i))

            bin_fthr = format(fthr, "03b")
            plt.xlabel("$\\tau_{p}$")
            plt.ylabel("rms [$DAC_{thr}$ code]")
            plt.title("Noise standard deviation with fthr = {}".format(bin_fthr))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.9, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.2, 0.5), borderaxespad=0, frameon=True, ncol=1)
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.tick_params(axis='x', which='minor', direction='out', bottom=True, length=3)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'B_fthr{}_TAU.svg'.format(fthr), format='svg', bbox_inches="tight")
            plt.close()

    return dac, a, b


"""=== DISPERSION MINIMISATION ===
Finding fthr which optimises threshold dispersion for each channel
"""


def dispersion_minimization(enablePlot, out_path_current, n_ch, n_tau, n_fthr, dac, a, b):
    cm = plt.get_cmap('hsv')
    colours = cm(np.linspace(0, 1.0, n_ch))

    # Output files
    file_dispersion = out_path_current + 'dispersion_minimization.dat'
    file_default = out_path_current + 'defaults_custom.defs'

    """=== PROCESSING ==="""
    print('Processing input files...')
    # with open(out_path_current + 'ThresholdScan_fthrOpt.dat','w') as f:
    # Outliers removal
    a_no_outliers = []
    b_outliers_inf = np.zeros((n_tau, n_ch * n_fthr))
    for i in range(n_ch):
        print('CH {}:'.format(i))
        a_caps = (plt.boxplot(a[:, i * n_fthr:i * n_fthr + n_fthr].flatten(), whis=[5, 95]))['caps']
        plt.close()
        b_caps = (plt.boxplot(b[:, i * n_fthr:i * n_fthr + n_fthr].flatten(), whis=[5, 95]))['caps']
        plt.close()
        min_a, max_a = a_caps[0].get_ydata()[0], a_caps[1].get_ydata()[0]
        min_b, max_b = b_caps[0].get_ydata()[0], b_caps[1].get_ydata()[0]
        print('  Values outside of interval [{},{}] for a are considered outliers.'.format(min_a, max_a))
        print('  Values outside of interval [{},{}] for b are considered outliers.\n'.format(min_b, max_b))

        for fthr in range(n_fthr):
            a_temp = np.copy(a[:, i * n_fthr + fthr])
            b_temp = np.copy(b[:, i * n_fthr + fthr])
            b_temp[np.where((b_temp <= min_b) | (b_temp >= max_b))] = math.inf
            a_no_outliers.append(a_temp[np.where((a_temp >= min_a) & (a_temp <= max_a))])
            b_outliers_inf[:, i * n_fthr + fthr] = b_temp

    print('Before outliers removal there were {} elements to consider for minimization.'.format(n_ch * n_tau * n_fthr))
    a_no_outliers_flat = np.array([item for a_no_outliers_fthr in a_no_outliers for item in a_no_outliers_fthr])
    print('After outliers removal there are {} elements to consider for minimization.'.format(len(a_no_outliers_flat)))
    a_mu_chip = np.mean(a_no_outliers_flat)

    vec_a_no_outliers_mu = np.zeros((n_fthr, n_ch))
    # Needs to be a loop cause each row has a different number of cols
    for i in range(n_ch):
        for fthr in range(n_fthr):
            vec_a_no_outliers_mu[fthr, i] = np.mean(a_no_outliers[i * n_fthr + fthr])

    vec_opt_fthr = np.zeros((n_ch), dtype=int)
    for i in range(n_ch):
        errors = np.abs(np.array(vec_a_no_outliers_mu[:, i]) - a_mu_chip)
        opt_fthr = np.argmin(errors)
        min_error = errors[opt_fthr]
        vec_opt_fthr[i] = opt_fthr
        # f.write('CH #{} optimized fthr={} with error={} wrt a of the chip (outliers excluded)'.format(i, opt_fthr, min_error))
        print('CH #{} optimized fthr={} with error={} wrt a of the chip (outliers excluded)'.format(i, opt_fthr, min_error))

    # Writes default.defs
    write_default(file_default, vec_opt_fthr, n_ch)

    # Writes output
    columns = ['fthr\\tp'] + [str(j) for j in range(n_tau)]
    data = np.zeros((2, 1 + n_tau))
    # For every channel @ fthr=000
    fthr = 0
    data[0, 0] = format(fthr, "03b")
    for j in range(n_tau):
        data[0, 1 + j] = np.std(a_no_outliers[fthr * n_tau::n_fthr][j])
    # For every channel tuned to optimised fthr
    for j in range(n_tau):
        for fthr in vec_opt_fthr:
            data[1, 1 + j] = np.std(a_no_outliers[fthr::n_fthr][j])

    np.savetxt(file_dispersion, data, delimiter='\t', comments='#{}'.format(','.join(columns)))

    if (enablePlot):
        """=== PLOTTING ==="""
        print('Plotting...')
        # Plot for all channel with optimisation
        print('Plotting optimized')
        for j in range(n_tau):
            print('  Tau {}:'.format(j))
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)
            for i in range(n_ch):
                x = np.linspace(100, 300, 500)
                bin_fthr = format(vec_opt_fthr[i], "03b")
                ax.plot(x, 100 * f_VTP(x, a[j, i * n_fthr + vec_opt_fthr[i]], b[j, i * n_fthr + vec_opt_fthr[i]]), c=colours[i], label="$CH \#{}, fthr={}$".format(i, bin_fthr))

            plt.xlabel("Discriminator threshold [$DAC_{thr}$ code]")
            plt.ylabel("Efficiency [%]")
            plt.title("Threshold scan with optimized fthr, for $\\tau_{{ {} }}$ and $DAC_{{inj}}$ code = {}".format(j, dac[0]))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.25, 0.5), borderaxespad=0, frameon=True, ncol=1)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'ThresholdScan_fthr_optimized_tau{}.svg'.format(j), format='svg', bbox_inches="tight")
            plt.close()

        # Plot for all channel without optimisation
        print('Plotting not optimized')
        bin_fthr = format(0, "03b")
        for j in range(n_tau):
            print('  Tau {}:'.format(j))
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)
            for i in range(n_ch):
                x = np.linspace(100, 300, 500)
                idx = i * n_fthr * n_tau + 0 * n_tau + j
                ax.plot(x, 100 * f_VTP(x, a[j, i * n_fthr + 0], b[j, i * n_fthr + 0]), c=colours[i], label="$CH \#{}$".format(i))

            plt.xlabel("Discriminator threshold [$DAC_{thr}$ code]")
            plt.ylabel("Efficiency [%]")
            plt.title("Threshold scan without optimization, fthr = {}, for $\\tau_{{ {} }}$ and $DAC_{{inj}}$ code = {}".format(bin_fthr, j, dac[0]))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.95, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.15, 0.5), borderaxespad=0, frameon=True, ncol=1)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current + 'ThresholdScan_fthr_not_optimized_tau{}.svg'.format(j), format='svg', bbox_inches="tight")
            plt.close()


"""=== FITTING FTHR ===
Same as in fitting TAU, but varying fthr
"""


def threshold_fitting_fthr(isProcessed, enablePlot, input_path, out_path_current, n_ch, n_tau, n_fthr, plt_fit_xmax, plt_fit_xmin):
    cm = plt.get_cmap('hsv')
    colours = cm(np.linspace(0, 1.0, n_ch))

    # Processed files
    sep = os.path.sep
    processed_path = out_path_current + 'Processed' + sep
    if (not(os.path.exists(processed_path))):
        os.makedirs(processed_path)
    file_a_FTHR = processed_path + 'a_fthr.dat'
    file_b_FTHR = processed_path + 'b_fthr.dat'

    a = np.zeros((n_tau, n_ch * n_fthr), dtype=float)
    b = np.zeros((n_tau, n_ch * n_fthr), dtype=float)
    dac = []

    if (isProcessed):
        a = np.loadtxt(file_a_FTHR, delimiter='\t', dtype=float)
        b = np.loadtxt(file_b_FTHR, delimiter='\t', dtype=float)
        try:
            fname = input_path + "ThresholdScan_fthr0_tau0.dat"
            dac = np.loadtxt(fname, dtype='int', comments='#', usecols=(1), unpack=True)
        except OSError:
            pass
    else:
        """=== PROCESSING ==="""
        print('Processing input files...')
        # with open(out_path_current + 'ThresholdScan.dat','w') as f:
        for i in range(n_ch):
            # f.write('CH #{}\n'.format(i))
            print('CH #{}'.format(i))
            for j in range(n_tau):
                plt.figure(figsize=(10, 6))
                ax = plt.subplot(111)
                #plt.xlim(np.floor(a_min) -5, np.floor(a_max) + 5)

                #f.write('  tau: {}\n'.format(j))
                print('  tau: {}'.format(j))
                for fthr in range(n_fthr):
                    #f.write('    fthr: {}\n'.format(fthr))
                    print('    fthr: {}'.format(fthr))
                    ok = 0
                    # Getting data
                    try:
                        fname = input_path + "ThresholdScan_fthr" + str(fthr) + "_tau" + str(j) + ".dat"
                        threshold, dac, events, triggered, ch = np.loadtxt(fname, dtype='int', comments='#', usecols=(0, 1, 2, 3, 4), unpack=True)
                    except OSError:
                        pass
                    thr = np.unique(threshold)
                    ev = np.unique(events)
                    y = np.zeros(len(thr))
                    ok = 1

                    # Plotting data
                    for k in range(len(thr)):
                        idx = (ch == i).nonzero()[0]
                        jdx = (threshold[idx] == thr[k]).nonzero()[0]
                        y[k] = triggered[idx[jdx]] / ev

                    # fitting function
                    b_init = np.random.normal(20, 1.5)
                    a_init = np.random.normal(250, 5)
                    x = thr
                    param = np.array([a_init, b_init])

                    try:
                        popt, pcov = optimize.curve_fit(f_VTP, x, y, param)  # bounds = (lower_bounds, upper_bounds))
                    except RuntimeError:
                        #f.write('    Not fitted for tau{}\n'.format(j))
                        print('    Not fitted for tau{}'.format(j))

                    # print(popt)
                    a[j, i * n_fthr + fthr] = popt[0]
                    b[j, i * n_fthr + fthr] = popt[1]

                    # Finding the threshold
                    #f.write('    Threshold for tau{}={:.2f}, b={:.2f}\n'.format(j, a_fthr[j],b_fthr[j]))
                    print('    Threshold for tau{}={:.2f}, b={:.2f}'.format(j, popt[0], popt[1]))

            #f.write('  Threshold for channel {}={:.2f}, sigma={:.2f}, dispersion={:.2f}\n'.format(i, a_mu, a_sigma, 3*2*a_sigma))
            a_mu, a_sigma = np.mean(a[:, i * n_fthr:i * n_fthr + n_fthr]), np.std(a[:, i * n_fthr:i * n_fthr + n_fthr])
            print('  Threshold for CH {}:'.format(i),
                  'mu={:.2f}'.format(a_mu),
                  ', sigma={:.2f}'.format(a_sigma),
                  ', dispersion:{:.2f}'.format(3 * 2 * a_sigma))

    np.savetxt(file_a_FTHR, a, delimiter='\t')
    np.savetxt(file_b_FTHR, b, delimiter='\t')

    if (enablePlot):
        """=== PLOTTING ==="""
        print('Plotting...')
        # with open(out_path_current + 'ThresholdScan.dat','w') as f:
        for i in range(n_ch):
            # f.write('CH #{}\n'.format(i))
            print('CH #{}'.format(i))
            for j in range(n_tau):
                plt.figure(figsize=(10, 6))
                ax = plt.subplot(111)
                #plt.xlim(np.floor(a_min) -5, np.floor(a_max) + 5)

                #f.write('  tau: {}\n'.format(j))
                print('  tau: {}'.format(j))
                for fthr in range(n_fthr):
                    #f.write('    fthr: {}\n'.format(fthr))
                    print('    fthr: {}'.format(fthr))

                    x = np.linspace(100, 300, 500)
                    bin_fthr = format(fthr, "03b")
                    ax.plot(x, 100 * f_VTP(x, a[j, i * n_fthr + fthr], b[j, i * n_fthr + fthr]), label="$CH \#{}, fthr={}$".format(i, bin_fthr))

                a_mu, a_sigma = np.mean(a[:, i * n_fthr:i * n_fthr + n_fthr]), np.std(a[:, i * n_fthr:i * n_fthr + n_fthr])
                ax.errorbar(a_mu, 50, xerr=np.multiply(3, a_sigma), xlolims=True, xuplims=True, c='b')
                ax.plot(a_mu, 50, c='b', marker='o', linewidth=2)

                plt.xlim(plt_fit_xmin, plt_fit_xmax)
                plt.xlabel("Discriminator threshold [$DAC_{thr}$ code]")
                plt.ylabel("Efficiency [%]")
                plt.title("Threshold of channel #{} with $DAC_{{inj}}$ = {} and $\\tau_{{p}} = {}$".format(i, dac[0], j))
                chartBox = ax.get_position()
                ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.75, chartBox.height])
                ax.legend(loc=7, bbox_to_anchor=(1.475, 0.5), borderaxespad=0, frameon=True, ncol=1)
                #plt.legend(loc = "lower right")
                plt.grid(True)
                # plt.show()
                plt.savefig(out_path_current + 'ThresholdScan_tau{}_ch{}.svg'.format(j, i), format='svg', bbox_inches="tight")
                plt.close()
