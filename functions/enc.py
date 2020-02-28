import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as plt_colours
import matplotlib.ticker as tickers
from matplotlib import cm

"""=== EQUIVALENT NOISE CHARGE ==="""
def enc(enablePlot, out_path_current, n_ch, n_tau, n_fthr, rms_pedestal, m_high_gain_lin, m_high_gain_poly_wide_range, poly_degree):
    out_string = ['_lin', '_poly']
    switch = [m_high_gain_lin, m_high_gain_poly_wide_range]
    title = ['linear', 'polynomial (deg={})'.format(poly_degree)]

    columns = ['tp\ch'] + [str(i) for i in range(n_ch)]
    data = np.zeros((n_tau, 1+n_ch))
    data[0:n_tau, 0] = [str(j) for j in range(n_tau)]

    #ENC
    num_bits_DAC = 16
    num_bits_ADC = 10
    C_inj = 1180e-15 #[fF]

    # First method:
    #   rms is the variance of the pedestal
    #   rms[ADC code] / gain [ADC code / DAC_inj code] -> [DAC_inj code]
    conv = 31.25 #[V / DAC_inj code]
    kn = 1.6e-19 / 3.6e-3 #[C / keV]
    gamma = (conv * C_inj * 10**(-6)) / kn #[kev / DAC_inj code]

    for s in range(len(switch)):
        """
        s:
            0 -> linear fit
            1 -> polynomial fit
        """


        for j in range(n_tau):
            for i in range(n_ch):
                if(switch[s][j, i] < 0.5):
                    switch[s][j, i] = np.nan


        enc_first_method = np.multiply(2.35 * gamma, np.divide(rms_pedestal, switch[s]))

        for i in range(n_ch):
            print('Channel #{}'.format(i))
            for j in range(n_tau):
                print('  tau_{}: {:.3f}'.format(j, enc_first_method[j, i]))
                data[j, i + 1] = enc_first_method[j, i]

        np.savetxt(out_path_current + 'enc_first_method' + out_string[s] + '.dat', data, fmt='%1.3f', delimiter='\t', comments='#{}'.format(','.join(columns)))

        if (enablePlot):
            # normal plot
            plt.figure(figsize=(15,10))
            ax = plt.subplot(111)
            x = np.linspace(0,31,32)

            for j in range(n_tau):
                ax.plot(x, enc_first_method[j], marker='*', markersize=5, linewidth=1, label='$\\tau_{}$'.format(j))

            xnew = np.linspace(0,32,1000)
            ax.plot(xnew, np.array([4 for i in range(len(xnew))]), linestyle='--', color='r', linewidth=1, label='_')
            plt.xlabel('Channels')
            plt.ylabel('FWHM ENC [keV]')
            plt.title('FWHM ENC [keV] (Pedestal rms / Gain with ' +  title[s] + ' fit slope)')
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.95, chartBox.height])
            ax.legend(loc=7, bbox_to_anchor=(1.08, 0.5), borderaxespad=0, frameon=True, ncol=1)
            ax.xaxis.set_minor_locator(tickers.AutoMinorLocator())
            ax.tick_params(axis='x',which='minor',direction='out',bottom=True, length=3)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            #plt.show()
            plt.savefig(out_path_current + 'enc_first_method' + out_string[s] + '.svg', format='svg', bbox_inches="tight")
            plt.close()

            # Barplot
            plt.figure(figsize=(20,10))
            ax_bar = plt.subplot(111)
            width = 0.09
            pos = [-7*width/2, -5*width/2, -3*width/2, -width/2, width/2, 3*width/2, 5*width/2, 7*width/2]

            for j in range(n_tau):
                ax_bar.bar(x + pos[j], enc_first_method[j], width=width, label='$\\tau_{}$'.format(j))

            plt.xlabel("Channels")
            plt.ylabel("FWHM ENC [keV]")
            plt.title('FWHM ENC [keV] (Pedestal rms / Gain with ' +  title[s] + ' fit slope)')
            chartBox = ax_bar.get_position()
            ax_bar.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.95, chartBox.height])
            ax_bar.legend(loc=7, bbox_to_anchor=(1.06, 0.5), borderaxespad=0, frameon=True, ncol=1)
            ax_bar.xaxis.set_minor_locator(tickers.AutoMinorLocator())
            ax_bar.tick_params(axis='x',which='minor',direction='out',bottom=True, length=3)
            #plt.legend(loc = "lower right")
            plt.grid(True)
            #plt.show()
            plt.savefig(out_path_current + 'enc_first_method' + out_string[s] + '_barplot.svg',format='svg', bbox_inches="tight")
            plt.close()

            # Heatmap
            plt.figure(figsize=(10,10))
            plt.set_cmap('jet')
            ax_heat = plt.subplot(111)

            current_cmap = cm.get_cmap()
            current_cmap.set_bad(color='purple')
            pos = ax_heat.imshow(enc_first_method.transpose())
            #ax_heat.set_xticks(np.arange(n_tau))
            #ax_heat.set_yticks(np.arange(n_ch))
            #ax_heat.set_xticklabels(np.linspace(0,7,8))
            #ax_heat.set_yticklabels(np.linspace(0,31,32))
            ax_heat.set_aspect('auto')
            for i in range(n_ch):
                for j in range(n_tau):
                    text = ax_heat.text(j, i, "%.2f" %  enc_first_method[j, i], ha='center', va='center')

            ax_heat.figure.colorbar(pos)
            plt.ylabel("Channels")
            plt.xlabel("$\\tau [\\mu s]$")
            plt.title('FWHM ENC [keV] (Pedestal rms / Gain with ' +  title[s] + ' fit slope)')
            plt.savefig(out_path_current + 'enc_first_method' + out_string[s] + '_heatmap.svg',format='svg', bbox_inches="tight")
            plt.close()

    # Second method: TODO
    #   rms is standard deviation of a [DAC_thr code]
    #   rms[DAC_thr code] / gain [ADC code / DAC_inj code]
    #gamma = 1
    #enc_second_method = np.multiply(2.35 *gamma, rms_b / gain_high)

    # Third method:
    #   rms is TODO
    #   rms[] / gain [ADC code / DAC_inj code]
    #gamma = 1
    #enc_third_method = np.multiply(2.35 * gamma, rms_ / gain_high)
