# Unit Match - Python

## Run UMPy

To run UMPy standard spike sorting data is needed; channel positions and 2 extracted raw waveforms for each unit. This can be calculated externally [BombCell](https://github.com/Julie-Fabre/bombcell) or using *extract_raw_demo.ipynb* to extract these waveforms from compressed data (.cbin and .ch) or raw data. There is also a [Spike Interface](https://spikeinterface.readthedocs.io/en/latest/) integrated notebook *UMPy_spike_interface_demo.ipynb* which uses spike interface to get this data.
Be careful not to mix and match the different ways of extracting raw waveforms, as there are difference between the methods.

There are to example notebooks for running UMPy *UMPy_example.ipynb* and *UMPy_example_detailed.ipynb*. These notebooks will guide you through running Unit Match all you need to supply is paths to the data. *UMPy_example.ipynb* is recommended to use first as it is simpler, however *UMPy_example_detailed.ipynb* may be useful in unique cases as it is more modular.

The GUI is a optional step to curated and investigate the information Unit Match has calculated; for efficient usage of the GUI please look at *GUI_Reference_Guide.md* in the Demo Notebooks folder.

## Dependencies

This version relies on many core python packages: [numpy](https://numpy.org/), [scipy](https://scipy.org/), [JobLib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/), [tkinter](https://docs.python.org/3/library/tkinter.html) and [matplotlib](https://matplotlib.org/). All of these libaries come with a [Anaconda](https://www.anaconda.com/download/) version of python. 
For extracting raw data, the library [mtscomp](https://github.com/int-brain-lab/mtscomp) is needed, and can be installed by `pip install mtscomp`.


## Installation 

After creating an python environment

```
conda create --name UnitMatch python==3.9 
conda activate UnitMatch
```

You can install UnitMatchPy with pip. It will automatically install all the dependencies.

```
pip install UnitMatchPy
```
