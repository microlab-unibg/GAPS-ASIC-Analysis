"""=== NOT IMPLEMENTED ==="""

"""## Equivalent Noise Charge"""
"""if (only == '' or only == 'enc'):
    #ENC
    #b is from the fit function
    #D_V is from the tf scan but is limited
    #gain_high is the gain of the transfer function low energy
    sup_voltage_DAC = 2.5
    sup_voltage_ADC = 1.8
    num_bits_DAC = 16
    num_bits_ADC = 10
    C_inj = 890e-15 #fF
    kn = 0.044 #[keV/C]
    gain_SH = 3.23

    #[muV]
    noise_V = b * (sup_voltage_DAC / (2 ** num_bits_DAC) * 10e6)

    #[keV/DACu]
    D_V = (sup_voltage_DAC / (2 ** num_bits_DAC))
    Q_inj = C_inj * kn * D_V

    #[muV/keV]
    gain_Vkev = gain_high * ((sup_voltage_ADC / (2 ** num_bits_ADC)) / Q_inj) * 10e6

    #[muV/keV]
    gain_shaper = gain_Vkev / gain_SH

    #[keV] MANTIENI LARGHEZZA MEZZA ALTEZZA
    enc = 2.35 * noise_V / gain_shaper

    print(enc)"""

""" Threshold scan with fixed tau and variable fthr
for i in range(n_ch):
    vec_th = []
    some_nans = 0

    for j in range(n_tau):
        ok = 0
        # Getting data
        try:
            fname = input_path + "ThresholdScan_tau4_ch" + str(i) + "_fthr" + str(j) + ".dat"
            threshold, dac, events, triggered = np.loadtxt(fname,dtype='int',comments='#',usecols=(0,1,2,3),unpack=True)
        except OSError:
            pass
        thr = np.unique(threshold)
        ev = np.unique(events)
        y = np.zeros(len(dac))
        ok = 1

        # Plotting data
        for k in range(len(dac)):
            y[k] = 100*triggered[k]/ev

        # Finding the threshold
        th_idxs = [i for i in range(len(y)) if y[i] >= 40 and y[i] <= 60]
        th = np.mean(dac[th_idxs])
        if (math.isnan(th)):
            some_nans = 1
        vec_th.append(th)
        print('Threshold for',j,'fthr:','{:.2f}'.format(th) if not(math.isnan(th)) else "no info")

        plt.plot(dac,y,label = "fthr = {}{}".format(j,""if not(math.isnan(th)) else "(no info)"))

    # Finding the threshold for the channel and the dispersion
    if (not(some_nans)):
        th_mu, th_var = stats.norm.fit(vec_th)
        print('Threshold for channel',i,':','{:.2f}'.format(th_mu),', var:','{:.2f}'.format(th_var),', dispersion:','{:.2f}'.format(3*2*th_var))

        plt.errorbar(th_mu,50,xerr=np.multiply(3,th_var),xlolims=True,xuplims=True,c='b')
        plt.plot(th_mu,50,c='b',marker='o',linewidth=2)

    if (ok):
        plt.xlabel("DAC code [AU]")
        plt.ylabel("efficiency [%]")
        plt.title("Threshold scan of channel #{} (threshold = {})".format(i,thr))
        plt.legend(loc = "lower right")
        plt.show()"""
