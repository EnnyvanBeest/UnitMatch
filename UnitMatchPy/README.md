# Unit Match - Python

## Run UMPy

To run UMPy, standard spike sorting data is needed; channel positions and 2 extracted raw waveforms for each unit. The raw waveforms can be extracted from compressed or raw data using [BombCell](https://github.com/Julie-Fabre/bombcell) or *ExtractRawDataDemo.ipynb*. There is also a [SpikeInterface](https://spikeinterface.readthedocs.io/en/latest/) integrated notebook *UMPy_spike_interface_demo.ipynb* which uses SpikeInterface to get this data.
Be careful not to mix and match the different ways of exrtacting raw waveforms, as there are differences between the methods.

There are two example notebooks for running UMPy: *UMPyExample.ipynb* and *UMPyExampleBrief.ipynb*. These notebooks will guide you through running UnitMatch, and all you need to supply is paths to the data. *UMPyExampleBried.ipynb* is recommended to use first as it is simpler. However, *UMPyExample.ipynb* is more modular and may thus be useful in specific cases.

The GUI is an optional step to curate and investigate the neuron pairs found by UnitMatch. For efficient usage of the GUI, please look at *GUI_Reference_Guide.md* in the Demo Notebooks folder.

## Dependencies

This version relies on many core python packages: [numpy](https://numpy.org/), [scipy](https://scipy.org/), [JobLib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/), [tkinter](https://docs.python.org/3/library/tkinter.html) and [matplotlib](https://matplotlib.org/). All of these libaries come with a [Anaconda](https://www.anaconda.com/download/) version of python. 
For extracting raw data, the library [mtscomp](https://github.com/int-brain-lab/mtscomp) is needed, and can be installed by `pip install mtscomp`.


## Installation 

After creating an python environment

```
conda create --name UnitMatch python 
conda activate UnitMatch
```

You can install UnitMatchPy with pip. It will automatically install all the dependencies.

```
pip install UnitMatchPy
```
