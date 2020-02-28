import numpy as np
from scipy import special

"""=== FUNCTIONS ==="""
def f_VTP(x, a, b):
    """Function to be fitted for threshold scan analysis
    Parameters
    ----------
    x : array
        Function input
    a : scalar, array
        Threshold voltage (to be fitted)
    b : scalar, array
        std dev of the threshold voltage (to be fitted)
    Returns
    -------
    Function evaluation
    """
    return (1/2 + 1/2 * special.erf( (x - a) / (np.sqrt(2) * b) ) )


def cubic_function(x, a, b, c, d):
    """Cubic function
    """
    return a * x**3 + b * x**2 + c * x + d


def write_data(file_name, columns, data):
    """Writes data in a table-formatted fashion in a csv file
    Parameters
    ----------
    file_name : str
        Name of the file to write to
    columns : list (n)
        List of strings which represent the first row of the table (metadata)
    data : Matrix (m x n)
        Matrix of data
    Returns
    -------
    None
    """
    np.savetxt(file_name, data, delimiter='\t', comments='#' + ','.join(columns))


def adjacent_values(vals, q1, q3):
    """Helper function for violinplot visualisation (courtesy of
    https://matplotlib.org/gallery/statistics/customized_violin.html#sphx-glr-gallery-statistics-customized-violin-py)
    """
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def write_default(file, vec_fthr, n_ch):
    """Writes default.defs file
    Parameters
    ----------
    file : str
        Path where to write defaults.defs
    vec_fthr : list
        Fine threshold for each channel
    Returns
    -------
    None
    """
    s = """[FPGA]
spi clock = 12
adc clock = 24
timeout = 500
event delay = 2000

[ASIC]
write adress = Broadcast
read adress = 1
mode = 000
bias = 0100
csa refs = 1100
shaper = 100
leakage mask = 00000000000000000000000000000000
trigger enable mask = 11111111111111111111111111111111
calibrazion mask = 11111111111111111111111111111111
threshold = 128
"""

    for i in range(n_ch):
        s = s + 'fine threshold ch#{} = {}\n'.format(i, format(int(vec_fthr[i]), "03b"))
    s = s + """

[Manual test]
inject = 0
count events = 0
dac = 100
dac max = 1000
step = 100
sweep dac = 0
hold delay = 48
sweep hold delay = 0
thr max = 200
thr min = 100
thr step = 1
sweep thr = 0
events = 10

[Automated test]
fast test = 1
enable pedestals = 1
pedestal events = 1000
enable waveform scan = 1
waveform scan events = 16
waveform scan dac = 1000
waveform scan max delay = 120
enable transfer function = 1
transfer function events = 16
range #1 dac min = 10
range #1 dac max = 100
range #1 step = 10
range #2 dac min = 200
range #2 dac max = 1000
range #2 step = 100
enable range #2 = 1
range #3 dac min = 1200
range #3 dac max = 2000
range #3 step = 200
enable range #3 = 1
range #4 dac min = 4000
range #4 dac max = 64000
range #4 step = 2000
enable range #4 = 1
enable threshold scan = 1
threshold scan events = 1000
threshold scan thr min = 200
threshold scan thr max = 255
threshold scan step = 1
enable self-trigger test = 1
self-trigger test events = 100
self-trigger test dac = 1000
self-trigger test threshold = 200
self-trigger test shaper = 4
enable adc test = 1
adc test events = 16
adc test dac min = 15000
adc test dac max = 42000
adc test step = 10
"""

    with open(file, 'w') as f:
        f.write(s)
