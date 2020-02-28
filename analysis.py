""" ANALYSIS
[Matteo Fratus, m.fratus(at)studenti.unibg.it]
[Paolo Lazzaroni, p.lazzaroni(at)studenti.unibg.it]

This script provides a little bit of insight into data that comes from ASICs under test.
The steps involved in the main computation are:
- Pedestal, where dispersion of a channel with Q_inj = 0 is measured and queried
- Waveform, which provides stats on shaper output waveform
- Transfer Function, in which high and low energy gains are computed
- Threshold, where threshold of channels are found as well as stats like dispersion,
    best fthr combination, etc
- ENC (to be implemented), evaluation of enc
Each section is detailed in the corresponding file under "functions" folder, which must be
in the same directory as this file is.

=== USAGE ===
A path is required when launching the script. This must be the path to a folder
which contains a subfolder named "data" with all needed .dat files. This script
will create a subfolder named "analysis" in the path (at the same level of "data"
subfolder) where output files will be stored.
The path must be relative to this script.
I.E.:
    Folder structure:
        - analysis.py
        - Measures
            > 00
            > 01
                * data
                    # *.dat
            > 02
    Relative path:
        Measures/01

If singleMeasures, explained under "PARAMETERS" in this introduction, is set to True, a
peaking time will be required and the analysis will be performed only wrt that peaking time.

=== PARAMETERS ===
only:
    This parameter has to be set to the step one wants to perform. Accepted values for this
    parameter are:
        '' (everything), pedestal, waveform, transfer, threshold, fitting_TAU,
        dispersion_min, fitting_FTHR, enc

isProcessed:
    If this parameter is set to True, the script will assume that input files have already
    been processed one time. If that is the case, while running, the script yielded some files
    (under folders named "Processed" in each analysis) which contains the results of the
    analysis and so it is to no use to compute again those metrics.

    This is useful for a ENC calculation workflow, as an example. One might need to analyse
    pedestal and transfer function and THEN decide to compute ENC. Since ENC needs pedestal and
    transfer function data, if those analysis have already been done, one might set isProcessed
    to True, so that the script will read processed data from folder "Processed" in Pedestal and
    Transfer function and skip those computation, getting straight to ENC.

    If this parameter is set to False, the script will process data assuming no prior analysis
    was performed.

singleMeasures:
    If this parameter is set to True, the script will ask for a single (for now) peaking time to
    be analysed, while other peaking times will be ignored. This parameter is supported for ENC analysis so far.

plotX (where X = {plotPedestal, Waveform, TransferFunction, ThresholdScanRaw, ThresholdScanFittedTau, ThresholdScanFittedFthr})
    These parameters will decide whether plots will be yielded for each analysis.
    If set to True, plots will be produced.
    If set to False, plots will not be produced.
"""

"""=== LIBRARIES ==="""
from functions.enc import enc
from functions.pedestal import pedestal
from functions.waveform import waveform
from functions.transfer import transfer
from functions.threshold import *
import os
import pathlib

# Global variables
n_ch = 32 # number of channels
n_tau = 8  # number of tau (peaking times)
n_fthr = 8  # number of fine thrimming settings
alpha = 1e-3  # confidence level for normality test

"""=== PARAMETERS ==="""
"""Which analysis is to be performed:
'' (everything), pedestal, waveform, transfer, threshold, fitting_TAU,
dispersion_min, fitting_FTHR, enc"""
only = 'enc'

"""Whether raw data are fetched from files and are to be processed or not
If True raw data files reading will be skipped as well as calculations, only
graphs will be updated. In some analysis, like 'waveform' and 'threshold' it is
not possible to skip computation, as one needs all data for plotting them.
If False everything will be computed as in a normal workflow."""
isProcessed = True

"""If this parameter is set to True, the script will expect to find some -
not every - files (for example wrt a single tau for each channel). The analysis
will be the same, input files will be searched in a different way (for now file name
will be checked, if the name std changes, check for the mask)"""
singleMeasures = True

"""Whether plots are plotted or not (disable for faster computation)"""
plotPedestal = True
plotWaveform = True
plotTransferFunction = True
plotThresholdScanRaw = True
plotThresholdScanFittedTau = True
plotThresholdScanFittedFthr = True

"""Test Parameters"""
# pedestal_events = 1000
# waveform_delay = (0, 119)
# transfer_dac_inj = (10, 64000)
# threshold_threshold = (200, 255)
#
# waveform_sample = waveform_delay[1] - waveform_delay[0]
# transfer_sample = transfer_dac_inj[1] - transfer_dac_inj[0]
# threshold_sample = threshold_threshold[1] - threshold_threshold[0]

"""=== MAIN COMPUTATION ==="""
# Getting where data are stored
print("Insert path to ASIC number folder (relative):")
sep = os.path.sep
path = str(pathlib.Path().absolute()) + sep + input()
if (not(os.path.isdir(path))):
    raise OSError("File not found")
    exit()

peaking_time = False
if (singleMeasures):
    print("Insert peaking time for single measure:")
    peaking_time = int(input())

# Data should be under a folder called 'data'
input_path = path + sep + 'data' + sep
if (not(path.endswith(sep))):
    path = path + sep
# Analysis outputs must be under a folder called 'analysis'
out_path = path + 'analysis' + sep
if (not(os.path.exists(out_path))):
    os.makedirs(out_path)
# Visual check
print(path, input_path, out_path)


""" === PEDESTAL ==="""
if (only == '' or only == 'pedestal' or only == 'enc'):
    out_path_current = out_path + 'Pedestal' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)

    rms_pedestal = pedestal(isProcessed, plotPedestal, singleMeasures, peaking_time, input_path, out_path_current, n_ch, n_tau, n_fthr)


""" === WAVEFORM === """
if (only == '' or only == 'waveform'):
    out_path_current = out_path + 'WaveformScan' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)

    waveform(plotWaveform, input_path, out_path_current, n_ch, n_tau, n_fthr)


"""=== TRANSFER FUNCTION === """
if (only == '' or only == 'transfer' or only == 'enc'):
    out_path_transfer_function = out_path + 'Transfer_Function' + sep
    out_path_current_CH = out_path_transfer_function + 'Transfer_Function_CH' + sep
    out_path_he = out_path_transfer_function + 'Low_Energy' + sep
    out_path_le = out_path_transfer_function + 'High_Energy' + sep
    if (not(os.path.exists(out_path_transfer_function))):
        os.makedirs(out_path_transfer_function)
    if (not(os.path.exists(out_path_current_CH))):
        os.makedirs(out_path_current_CH)
    if (not(os.path.exists(out_path_he))):
        os.makedirs(out_path_he)
    if (not(os.path.exists(out_path_le))):
        os.makedirs(out_path_le)

    out_path_current_TAU = out_path_transfer_function + 'Transfer_Function_TAU' + sep
    if (not(os.path.exists(out_path_current_TAU))):
        os.makedirs(out_path_current_TAU)

    m_high_gain_lin, m_high_gain_poly_wide_range, degree_fit = transfer(isProcessed, plotTransferFunction, singleMeasures, peaking_time, input_path, out_path_current_CH, out_path_he, out_path_le, out_path_current_TAU, n_ch, n_tau, n_fthr)


"""=== THRESHOLD SCAN === """
""" Raw analysis """
# False is for skipping this part
if (False and only == '' or only == 'threshold'):
    out_path_current = out_path + 'Threshold_Scan_TAU' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)

    threshold(plotThresholdScanRaw, input_path, out_path_current, n_ch, n_tau, n_fthr)


""" Min and max for pretty plotting """
if (False and only == '' or only == 'fitting_TAU' or only == 'dispersion_min' or only == 'fitting_FTHR'):
    plt_fit_xmax, plt_fit_xmin = find_min_max_threshold()


""" Fitting function for threshold scan by tau """
if (False and only == '' or only == 'fitting_TAU' or only == 'dispersion_min'):
    out_path_current = out_path + 'Threshold_Scan_Fitted_TAU' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)

    dac, a, b = threshold_fitting_tau(isProcessed, plotThresholdScanFittedTau, input_path, out_path_current, n_ch, n_tau, n_fthr, plt_fit_xmax, plt_fit_xmin)


""" Find set of fthrs which minimises dispersion """
if (False and only == '' or only == 'dispersion_min'):
    out_path_current = out_path + 'Threshold_Scan_Fthr_Optimization' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)

    # Needs some revision, does not plot the right graphs
    dispersion_minimization(True, out_path_current, n_ch, n_tau, n_fthr, dac, a, b)


""" Fitting function for threshold scan by fthr """
if (False and only == '' or only == 'fitting_FTHR'):
    out_path_current = out_path + 'Threshold_Scan_Fitted_FTHR' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)

    threshold_fitting_fthr(isProcessed, plotThresholdScanFittedFthr, input_path, out_path_current, n_ch, n_tau, n_fthr, plt_fit_xmax, plt_fit_xmin)

"""## Equivalent Noise Charge"""
if (only == '' or only == 'enc'):
    out_path_current = out_path + 'enc' + sep
    if (not(os.path.exists(out_path_current))):
        os.makedirs(out_path_current)
    enc(True, out_path_current, n_ch, n_tau, n_fthr, rms_pedestal, m_high_gain_lin, m_high_gain_poly_wide_range, degree_fit)
