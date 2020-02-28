# GAPS ASIC Analysis
This script provides a little bit of insight into data that comes from ASICs under test.
The steps involved in the main computation are:
- Pedestal, where dispersion of a channel with Q_inj = 0 is measured and queried
- Waveform, which provides stats on shaper output waveform
- Transfer Function, in which high and low energy gains are computed
- Threshold, where threshold of channels are found as well as stats like dispersion,
    best fthr combination, etc
- ENC, evaluation of ENC
Each section is detailed in the corresponding file under "functions" folder, which must be
in the same directory as this file is.

## Usage
A path is required when launching the script. This must be the path to a folder
which contains a subfolder named `data` with all needed `.dat` files. This script
will create a subfolder named `analysis` in the path (at the same level of `data`
subfolder) where output files will be stored.
The path must be relative to this script.
I.E.: Folder structure:
- analysis.py
- Measures
  - 00
  - 01
    - data
      - \*.dat
  - 02
Relative path: `Measures/01`

If `singleMeasures`, explained under "Parameters" in this introduction, is set to `True`, a
peaking time will be required and the analysis will be performed only wrt that peaking time.

## Parameters
`only`:
    This parameter has to be set to the step one wants to perform. Accepted values for this
    parameter are:
        '' (everything), pedestal, waveform, transfer, threshold, fitting_TAU,
        dispersion_min, fitting_FTHR, enc

`isProcessed`:
    If this parameter is set to `True`, the script will assume that input files have already
    been processed one time. If that is the case, while running, the script yielded some files
    (under folders named `Processed` in each analysis) which contains the results of the
    analysis and so it is to no use to compute again those metrics.

    This is useful for a ENC calculation workflow, as an example. One might need to analyse
    pedestal and transfer function and THEN decide to compute ENC. Since ENC needs pedestal and
    transfer function data, if those analysis have already been done, one might set `isProcessed`
    to `True`, so that the script will read processed data from folder `Processed` in Pedestal and
    Transfer function and skip those computation, getting straight to ENC.

    If this parameter is set to `False`, the script will process data assuming no prior analysis
    was performed.

`singleMeasures`:
    If this parameter is set to `True`, the script will ask for a single (for now) peaking time to
    be analysed, while other peaking times will be ignored.

`plotX` (where X = {plotPedestal, Waveform, TransferFunction, ThresholdScanRaw, ThresholdScanFittedTau, ThresholdScanFittedFthr})
    These parameters will decide whether plots will be yielded for each analysis.
    If set to `True`, plots will be produced.
    If set to `False`, plots will not be produced.

## Authors
- Matteo Fratus, m.fratus(at)studenti.unibg.it
- Paolo Lazzaroni, p.lazzaroni(at)studenti.unibg.it

## License
MIT License

Copyright (c) 2020 Matteo Fratus and Paolo Lazzaroni
