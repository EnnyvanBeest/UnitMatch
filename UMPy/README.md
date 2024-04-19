# Unit Match - Python

## Run UMPy

To run UMPy standard spike sorting data is needed; channel positions and 2 extractracted raw waveforms for each unit. This can be calculated exteranlly [BombCell](https://github.com/Julie-Fabre/bombcell) or using *ExtractRawDataDemo.ipynb* to extract these waveforms from compressed data (.cbin and .ch) or raw data. There is also a [Spike Interface](https://spikeinterface.readthedocs.io/en/latest/) intergrated notebook *UMPy_spike_interface_demo.ipynb* which uses spike interface to get this data.
Be careful not to mix and match the different ways of exrtacting raw waveforms, as there are difference between the methods.

There are to example notebooks for running UMPy *UMPyExample.ipynb* and *UMPyExampleBrief.ipynb*. These notebooks will guide you throuh running Unit Match all you need to supply is paths to the data. *UMPyExampleBried.ipynb* is reccomened to use first as it is simpler, however *UMPyExample.ipynb* may be useful in unique cases as it is more modular.

The GUI is a optional step to curated and investigate the infomation Unit Match has calculated; for efficent usage of the GUI please look at *GUI_Reference_Guide.md* in the Demo Notebooks folder.
## Dependencies
This version relies on many core python packages [numpy](https://numpy.org/), [scipy](https://scipy.org/), [JobLib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/), [tkinter](https://docs.python.org/3/library/tkinter.html) and [matplotlib](https://matplotlib.org/). All of these libaries come with a [Anaconda](https://www.anaconda.com/download/) version of python. 
For extracting raw data the library [mtscomp](https://github.com/int-brain-lab/mtscomp) is needed, this can be installed by `pip install mtscomp`.


## Installation 

After creating an python environment

```
conda create --name umpy python 
conda activate umpy
```

You can install UMPy with pip. It will automatically install all the dependency.

```
pip install -e /path_to_UMPy/UMPy
```
