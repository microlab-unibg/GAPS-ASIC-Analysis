import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

import os

"""=== TRANSFER FUNCTION ===
As the name states

### Measures
- Kink
- Gain in high and low frequency"""


def transfer(isProcessed, enablePlot, singleMeasures, peaking_time, input_path, out_path_current_CH, out_path_he, out_path_le, out_path_current_TAU, n_ch, n_tau, n_fthr):
    list_tau = [peaking_time] if singleMeasures else np.arange(n_tau)

    for j in list_tau:
        fname = input_path + "TransferFunction_fast_tau{}.dat".format(j)
        break

    try:
        daca = np.loadtxt(fname, comments='#', usecols=(1), unpack=True)
    except OSError:
        print('No Pedestal file was found. Please check that pedestal file name meets "TransferFunction_fast_tauX.dat" pattern.')
        exit(1)

    sep = os.path.sep
    degree_fit = 4

    colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

    # Processed files
    processed_path = out_path_current_CH + 'Processed' + sep
    if (not(os.path.exists(processed_path))):
        os.makedirs(processed_path)
    processed_path_plot = processed_path + 'SupportData' + sep
    if (not(os.path.exists(processed_path_plot))):
        os.makedirs(processed_path_plot)

    file_m_y = processed_path_plot + 'y.dat'
    file_high_gain_lin = processed_path + 'high_gain_lin.dat'
    file_high_gain_poly_wide_range = processed_path + 'high_gain_poly_wide_range.dat'
    file_low_gain_lin = processed_path + 'low_gain_lin.dat'
    file_lin_intercept_high = processed_path_plot + 'lin_intercept_high.dat'
    file_lin_intercept_low = processed_path_plot + 'lin_intercept_low.dat'
    file_poly_coeff_wide = processed_path_plot + 'poly_coeff_wide.dat'

    # Processed
    dac = np.unique(daca)
    m_y = np.zeros((n_tau, n_ch * len(dac)), dtype='float')
    m_high_gain_lin = np.zeros((n_tau, n_ch), dtype='float')
    m_high_gain_poly_wide_range = np.zeros((n_tau, n_ch), dtype='float')
    m_low_gain_lin = np.zeros((n_tau, n_ch), dtype='float')
    m_lin_intercept_high = np.zeros((n_tau, n_ch), dtype='float')
    m_lin_intercept_low = np.zeros((n_tau, n_ch), dtype='float')
    m_poly_coeff_wide = np.zeros((n_tau, n_ch * (degree_fit + 1)), dtype='float')

    if (isProcessed):
        try:
            m_y = np.loadtxt(file_m_y, dtype='float', delimiter='\t')
            m_high_gain_lin = np.loadtxt(file_high_gain_lin, dtype='float', delimiter='\t')
            m_high_gain_poly_wide_range = np.loadtxt(file_high_gain_poly_wide_range, dtype='float', delimiter='\t')
            m_low_gain_lin = np.loadtxt(file_low_gain_lin, dtype='float', delimiter='\t')
            m_lin_intercept_high = np.loadtxt(file_lin_intercept_high, dtype='float', delimiter='\t')
            m_lin_intercept_low = np.loadtxt(file_lin_intercept_low, dtype='float', delimiter='\t')
            m_poly_coeff_wide = np.loadtxt(file_poly_coeff_wide, dtype='float', delimiter='\t')
        except OSError:
            print('There are no processed files or some are missing for this analysis')
            exit(1)

    else:
        """=== PROCESSING ==="""

        # Data write tau for each channel
        data_gain_low = np.zeros((n_tau, 1 + n_ch))
        data_gain_high_lin = np.zeros((n_tau, 1 + n_ch))
        data_gain_high_poly = np.zeros((n_tau, 1 + n_ch))

        for i in range(n_ch):
            print('CH #{}'.format(i))
            for j in list_tau:
                print("  tau: {}".format(j))
                try:
                    fname = input_path + "TransferFunction_fast_tau" + str(j) + ".dat"
                    daca, typea, ch, val = np.loadtxt(fname, comments='#', usecols=(1, 2, 3, 4), unpack=True)
                    print('File TransferFunction_fast_tau{}.dat found!'.format(j))
                except OSError:
                    print('File TransferFunction_fast_tau{}.dat not found! Ignoring...'.format(j))
                    pass

                val[np.where((typea == 1) | (typea == 11))] = np.nan
                y = np.zeros(len(dac))
                for k in range(len(dac)):
                    idx = (ch == i).nonzero()[0]
                    jdx = (daca[idx] == dac[k]).nonzero()[0]
                    y[k] = np.nanmean(val[idx[jdx]])
                    m_y[j, i * len(dac) + k] = y[k]

        for i in range(n_ch):
            print('#CH {}'.format(i))
            for j in list_tau:
                print('tau {}'.format(j))
                # Low energy
                # linear interpolation [50, 100]
                base_idx = i * len(dac)
                this_y = m_y[j, base_idx:base_idx + len(dac)]
                dac_le = dac[np.where((dac >= 50) & (dac <= 100))]
                y_le = this_y[np.where((dac >= 50) & (dac <= 100))]
                (high_gain, intercept, r_value, p_value, std_err) = stats.linregress(dac_le, y_le)
                m_lin_intercept_high[j, i] = intercept
                m_high_gain_lin[j, i] = high_gain
                print('    Low Enegry gain lin: {:.3f}'.format(high_gain))
                print('      r_value: {:.3f}'.format(r_value))

                # polynomial interpolation with weighted initial points and derivative extraction [10 - 500]
                dac_le = dac[np.where(dac <= 500)]
                y_le = this_y[np.where(dac <= 500)]
                try:
                    popt = np.polyfit(dac_le[1:], y_le[1:], deg=degree_fit)
                except RuntimeError:
                    print('    Not fitted for tau{}'.format(j))

                m_high_gain_poly_wide_range[j, i] = popt[degree_fit - 1]
                for d in range(degree_fit + 1):
                    m_poly_coeff_wide[j, i * (degree_fit + 1) + d] = popt[d]

                poly = np.poly1d(popt)

                print('    Low Enegry gain poly: {:.3f}'.format(popt[degree_fit - 1]))
                ss_res = np.sum((y_le - poly(dac_le))**2)
                ss_tot = np.sum((y_le - np.mean(y_le))**2)
                r_value = 1 - (ss_res / ss_tot)
                print('      r_value: {:.3f}'.format(r_value))

                # High energy
                # linear interpolation [20000, 60000]
                dac_he = dac[np.where((dac >= 20000) & (dac <= 60000))]
                y_he = this_y[np.where((dac >= 20000) & (dac <= 60000))]
                if y_he.any():
                    (low_gain, intercept, r_value, p_value, std_err) = stats.linregress(dac_he, y_he)
                    m_low_gain_lin[j, i] = low_gain
                    m_lin_intercept_low[j, i] = intercept
                    print('    High energy gain: {:.3f}'.format(low_gain))
                    print('      r_value: {:.3f}'.format(r_value))

        # Save processed data
        header = 'Peak times (rows) and Channels, DAC (cols)\n'
        header = header + '\t'.join(['ch{}-dac{}'.format(i, daci) for i, daci in zip(range(n_ch), range(len(dac)))])
        np.savetxt(file_m_y, m_y, delimiter='\t', header=header)
        header_tau_ch = 'Peak times (rows) and Channels (cols)\n'
        header_tau_ch = header_tau_ch + '\t'.join(['ch{}'.format(i) for i in range(n_ch)])
        np.savetxt(file_high_gain_lin, m_high_gain_lin, delimiter='\t', header=header_tau_ch)
        np.savetxt(file_high_gain_poly_wide_range, m_high_gain_poly_wide_range, delimiter='\t', header=header_tau_ch)
        np.savetxt(file_low_gain_lin, m_low_gain_lin, delimiter='\t', header=header_tau_ch)
        np.savetxt(file_lin_intercept_high, m_lin_intercept_high, delimiter='\t', header=header_tau_ch)
        np.savetxt(file_lin_intercept_low, m_lin_intercept_low, delimiter='\t', header=header_tau_ch)
        header = 'Peak times (rows) and Channels, poly coefficient (cols)\n'
        header = header + '\t'.join(['ch{}-coef{}'.format(i, polyi) for i, polyi in zip(range(n_ch), range(degree_fit))])
        np.savetxt(file_poly_coeff_wide, m_poly_coeff_wide, delimiter='\t', header=header)

        print('Chip')
        header = 'DAC (rows) and Channels (cols)\n'
        header = header + 'dac\t' + '\t'.join(['ch{}'.format(i, polyi) for i, polyi in zip(range(n_ch), range(degree_fit))])
        for j in list_tau:
            if (y_he.any()):
                gain_mu, gain_sigma = np.mean(m_low_gain_lin[j]), np.std(m_low_gain_lin[j])
                #print('  tau ','{}'.format(i), ' Low gain mean:','{:.3f}'.format(j, gain_mu),', Low gain sigma','{:.3f}'.format(gain_sigma))
                print('  tau {}, Low gain mean: {:.3f}, Low gain sigma: {:.3f}'.format(j, gain_mu, gain_sigma))
            gain_mu, gain_sigma = np.mean(m_high_gain_lin[j]), np.std(m_high_gain_lin[j])
            print('  tau {}, High gain lin mean: {:.3f}, High gain lin sigma: {:.3f}'.format(j, gain_mu, gain_sigma))
            gain_mu, gain_sigma = np.mean(m_high_gain_poly_wide_range[j]), np.std(m_high_gain_poly_wide_range[j])
            print('  tau {}, High gain poly mean: {:.3f}, High gain poly sigma: {:.3f}'.format(j, gain_mu, gain_sigma))

            # Save data for each tau, values of ADC for each channel
            file_tau_ch = out_path_current_TAU + 'Values_tf_allch_tau' + str(j) + '.dat'
            m_out = np.zeros((len(dac), n_ch + 1), dtype='float')
            m_out[:, 0] = dac.transpose()
            for i in range(n_ch):
                m_out[:, i + 1] = m_y[j, i * len(dac):i * len(dac) + len(dac)].transpose()
            np.savetxt(file_tau_ch, m_out, delimiter='\t', header=header_tau_ch)

    if (enablePlot):
        """=== PLOTTING ==="""
        print('Plotting')

        for i in range(n_ch):
            print('  #CH {}'.format(i))
            y = np.zeros((1, len(dac)), dtype='float')

            fig, ax = plt.subplots(figsize=(10, 6))
            fig3, ax3 = plt.subplots(figsize=(10, 6))

            out_path_single = out_path_he + 'Channel_' + str(i) + sep
            if (not(os.path.exists(out_path_single))):
                os.makedirs(out_path_single)

            for j in list_tau:
                print('    tau {}'.format(j))
                fig2, ax2 = plt.subplots(figsize=(12, 6))
                y = m_y[j][i * len(dac):i * len(dac) + len(dac)]

                popt = np.array(m_poly_coeff_wide[j][i * (degree_fit + 1): i * (degree_fit + 1) + degree_fit + 1])

                ax.set_ylim(0, 2250)
                ax.plot(dac, y, label='$\\tau_{}$'.format(j))

                #spl = interpolate.UnivariateSpline(thightdac_le, y_le, s=1000, k=1)
                xnew = np.linspace(0, 500, 1000, endpoint=True)
                ax2.set_ylim(100, 500)
                #ax2.plot(dac_le, y_le, 'o', xnew, spl(xnew), '--', label='$\\tau_{}, G_0 = {:.3f}$'.format(j, gain_high), c=colours[j])
                ax2.plot(dac[np.where(dac <= 500)], y[np.where(dac <= 500)], 'o', mfc='none', label='data', color='k')
                #p2, = ax2.plot(xnew, spl(xnew), '--', label='_', c=colours[j])
                ax2.plot(xnew, m_lin_intercept_high[j][i] + m_high_gain_lin[j][i] * xnew, '-', label='linear interpolation, $G_0$ = {:.3f}'.format(m_high_gain_lin[j][i]), color='b')

                #ax2.plot(xnew, cubic_function(xnew, popt[0], popt[1], popt[2], popt[3]), '--', label='cubic interpolation [50-200], $G_0$ = {:.3f}'.format(popt[2]), color='r')
                #ax2.plot(xnew, cubic_function(xnew, popt_2[0], popt_2[1], popt_2[2], popt_2[3]), '-.', label='cubic interpolation [50-500], $G_0$ = {:.3f}'.format(popt_2[2]), color='g')

                poly = np.poly1d(popt)
                ax2.plot(xnew, poly(xnew), '-.', label='power 4 interpolation [10-500], $G_0$ = {:.3f}'.format(popt[degree_fit - 1]), color='r')

                #ax2.plot(-10,-10, label='$\\tau_{}, G_0 = {:.3f}$'.format(j, gain_high), linestyle='--', marker='o', c=colours[j])
                plt.xlabel("Cal_Voltage [$DAC_{inj}$ code]")
                plt.ylabel("Channel_out [ADC code]")
                plt.title("Low energy gain for channel #{}, tau: {}".format(i, j))
                chartBox = ax2.get_position()
                ax2.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.7, chartBox.height])
                ax2.legend(loc=7, bbox_to_anchor=(1.55, 0.5), borderaxespad=0, frameon=True, ncol=1)
                #plt.legend(loc = "lower right")
                plt.grid(True)
                # plt.show()
                plt.savefig(out_path_single + 'TransferFunctionLowEnergy_ch' + str(i) + '_tau' + str(j) + '.svg', format='svg', bbox_inches="tight")
                plt.close()

                # High energy
                dac_he = dac[np.where((dac >= 20000) & (dac <= 60000))]
                y_he = y[np.where((dac >= 20000) & (dac <= 60000))]
                if (y_he.any()):
                    xnew = np.linspace(dac_he[0], dac_he[-1], num=len(y_he), endpoint=True)
                    p1, = ax3.plot(dac_he, y_he, 'o', label='_', c=colours[j])
                    p2, = ax3.plot(xnew, m_lin_intercept_low[j][i] + m_low_gain_lin[j][i] * xnew, '--', label='_', c=colours[j])
                    ax3.plot(0, 0, label='$\\tau_{}, G_0 = {:.3f}$'.format(j, m_low_gain_lin[j][i]), linestyle='--', marker='o', c=colours[j])
                    # l = ax3.legend([(p1, p2)], ['$\\tau_{}, G_0 = {:.3f}$'.format(j, gain_low)], numpoints=1,
                    #                handler_map={tuple: legend_handler.HandlerTuple(ndivide=None)})
                    #plt.legend(['data', 'knots'], loc='best')
                    # plt.show()

            if (y_he.any()):
                plt.xlabel("Cal_Voltage [$DAC_{inj}$ code]")
                plt.ylabel("Channel_out [ADC code]")
                plt.title("High energy gain for channel #{}".format(i))
                chartBox = ax3.get_position()
                ax3.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.9, chartBox.height])
                ax3.legend(loc=7, bbox_to_anchor=(1.23, 0.5), borderaxespad=0, frameon=True, ncol=1)
                plt.ylim(1000, 2000)
                plt.xlim(19000, 61000)
                #plt.legend(loc = "lower right")
                plt.grid(True)
                # plt.show()
                plt.savefig(out_path_le + 'TransferFunctionHighEnergy_ch' + str(i) + '.svg', format='svg', bbox_inches="tight")
            plt.close()

            plt.xlabel("Cal_Voltage [$DAC_{inj}$ code]")
            plt.ylabel("Channel_out [ADC code]")
            plt.title("Transfer function of channel #{}".format(i))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.95, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.15, 0.5), borderaxespad=0, frameon=True, ncol=1)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current_CH + 'TransferFunction_ch' + str(i) + '.svg', format='svg', bbox_inches="tight")
            plt.close()

        # Plot channel for each tau
        for j in list_tau:
            plt.figure(figsize=(12, 11))
            ax = plt.subplot(111)
            ax.set_ylim(0, 2250)
            for i in range(n_ch):
                ax.plot(dac, m_y[j][i * len(dac): i * len(dac) + len(dac)], label='CH #{}'.format(i))

            plt.xlabel("Cal_Voltage [$DAC_{inj}$ code]")
            plt.ylabel("Channel_out [ADC code]")
            plt.title("Transfer function of $\\tau_{}$".format(j))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.95, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.15, 0.5), borderaxespad=0, frameon=True, ncol=1)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            # plt.show()
            plt.savefig(out_path_current_TAU + 'TransferFunction_tau' + str(j) + '_allch.svg', format='svg', bbox_inches="tight")
            plt.close()

    return m_high_gain_lin, m_high_gain_poly_wide_range, degree_fit
