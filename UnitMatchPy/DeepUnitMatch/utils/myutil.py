import numpy as np  
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import pandas as pd


def get_default_param(param = None):
    """
    Create param, a dictionary with the default parameters.
    If a dictionary is given, it will add values to it without overwriting existing values.
    Do not need to give a dictionary.
    """
    tmp = {'nTime' : 82, 'nChannels' : 384, 'ChannelRadius' : 110,
           'RnChannels' : 30, 'RnTime' : 60, 
        }
    # if no dictionary is given just returns the default parameters
    if param == None:
        out = tmp
    else:    
        # Add default parameters to param dictionary, does not overwrite pre existing param values
        out = tmp | param
    if out['RnChannels'] %2 !=0:
        print('RnChannels is not even, please check')
    return out

def get_threshold(mt:pd.DataFrame, metric:str="DNNSim", vis:bool=True, MAP=False):
    mt = mt.loc[(mt["RecSes1"]==mt["RecSes2"]), :]              # Only use within day rows to compute threshold
    n = np.sqrt(len(mt))

    # On-diagonal means same neuron. Off-diagonal means different neurons.
    on_diag = mt.loc[(mt["ID1"]==mt["ID2"]), [metric]]
    off_diag = mt.loc[(mt["ID1"]!=mt["ID2"]), [metric]]

    if len(on_diag) < 2 or len(off_diag) < 2:
        # can't do KDE in this case, or define a threshold properly.
        return mt[metric].max()

    # sanity check that the categories are being loaded correctly
    assert sum(on_diag.index.isin(off_diag.index)) == 0

    # Kernel density estimation (distributions are more useful than histograms)
    kde_on = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(on_diag.values)
    kde_off = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(off_diag.values)
    x = np.linspace(min(off_diag[metric]), max(on_diag[metric]), 1000).reshape(-1, 1)
    y_on = np.exp(kde_on.score_samples(x))
    y_off = np.exp(kde_off.score_samples(x))

    # Find the threshold where the distributions intersect
    if MAP:
        thresh = np.argwhere(np.diff(np.sign(y_off*(n-1) - y_on)))
    else:
        thresh=np.argwhere(np.diff(np.sign(y_off - y_on)))
    if len(thresh) == 0:
        thresh = len(x) - 1
    elif len(thresh) > 1:
        thresh = thresh[-1]

    if vis:
        # visualise the results
        print(f"Threshold: {metric} = {x[thresh].item()}")
        plt.hist(on_diag[metric], bins = 500, alpha = 0.5, density=True, label="On diagonal")
        plt.hist(off_diag[metric], bins = 500, alpha = 0.5, density=True, label="Off diagonal")
        plt.plot(x, y_on, label="On diagonal")
        plt.plot(x, y_off, label="Off diagonal")
        plt.axvline(x = x[thresh].item())
        
        plt.grid()
        plt.legend()
        plt.xlabel(metric)
        plt.title(f"Normalised histograms for on- and off-diagonal {metric} values in the same recording")
        plt.show()
    return x[thresh].item()

# def read_pos(path, skip_removed_units=False):
#     files = os.listdir(path)
#     x = []
#     y = []
#     filenames = []

#     for file in files:
#         if skip_removed_units:
#             if '+' in file or '#' in file:
#                 continue
#         fp = os.path.join(path, file)
#         with h5py.File(fp, 'r') as f:
#             # waveform = f['waveform'][()] 
#             MaxSitepos = f['MaxSitepos'][()]
#         x.append(MaxSitepos[0])
#         y.append(MaxSitepos[1])
#         filenames.append(get_unit_id(file))
#     output = {"file":filenames, "x":x, "y":y}

#     return output
