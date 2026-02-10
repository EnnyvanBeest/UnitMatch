# UnitMatchPy (UMPy) and DeepUnitMatch (DUM)

This repository contains:
- **UnitMatchPy (UMPy)**: automatic matching of neurons across sessions (Python package).
- **DeepUnitMatch (DUM)**: DeepUnitMatch pipeline and demos (see preprint below).

## References

- **UnitMatch (UMPy):** https://www.nature.com/articles/s41592-024-02440-1
- **DeepUnitMatch (DUM) preprint:** https://www.biorxiv.org/content/10.64898/2026.01.30.702777v1

## Versions

- **UnitMatchPy version:** `3.2.9` (from `pyproject.toml`)
- **DeepUnitMatch:** code lives under `DeepUnitMatch/` (see **DeepUnitMatch (DUM)** section)

## Installation

`pip install` (including `pip install -e .`) installs into whatever Python environment your `pip` points to (system Python, a conda env, or a virtualenv). `pip` does **not** create or name environments.

We recommend using Anaconda/Miniconda (conda) to create an isolated environment first:

```bash
# Create a new environment (pick any name you like; example: UMPy)
conda create -n UMPy python=3.11 pip #(press y when prompted)
conda activate UMPy
```

Then install using pip (options below).

### Option A: Install the released package (PyPI)

```bash
pip install UnitMatchPy
```

Optional extras (heavier dependencies used by some notebooks and integration with SpikeInterface):

```bash
pip install "UnitMatchPy[full,notebooks]"
```

### Option B: Install a local, editable copy (for development / modified code)

First, open a terminal and navigate to this folder (the one containing `pyproject.toml`). The `pip install -e` command must be run from here:

```bash
# Windows (PowerShell)
cd $HOME\Documents\GitHub\DeepUnitMatch\UnitMatchPy

# macOS / Linux
cd ~/Documents/GitHub/DeepUnitMatch/UnitMatchPy
```

```bash
pip install -e .
```

Optional extras (e.g. if you'd like to run the notebooks or integrate with SpikeInterface):

```bash
pip install -e ".[full,notebooks]"
```

## Demo notebooks

All demo notebooks are in `Demo Notebooks/`.

### Run UnitMatchPy (UMPy)

To run UMPy, standard spike sorting data is needed (channel positions and extracted raw waveforms for each unit). Waveforms can be extracted externally (e.g. [BombCell](https://github.com/Julie-Fabre/bombcell)) or using the demo notebooks:
- `Demo Notebooks/extract_raw_data_demo.ipynb` (compressed `.cbin`/`.ch` or raw)
- `Demo Notebooks/extract_raw_data_demo_open_ephys.ipynb` (Open Ephys)
- `Demo Notebooks/UMPy_spike_interface_demo.ipynb` ([SpikeInterface](https://spikeinterface.readthedocs.io/en/latest/) workflow)

Example notebooks:
- `Demo Notebooks/UMPy_example.ipynb` (recommended starting point)
- `Demo Notebooks/UMPy_example_detailed.ipynb` (more modular / advanced)

The GUI is an optional step to curate and explore UnitMatch outputs; see `Demo Notebooks/GUI_Reference_Guide.md` for usage tips and shortcuts.

## DeepUnitMatch (DUM)

To try the DeepUnitMatch version, start with:
- `Demo Notebooks/DeepUnitMatch.ipynb`

If you want to train / fine-tune a model on your own data, see:
- `Demo Notebooks/DUM_training.ipynb`

Important note: DeepUnitMatch current trained model is for Npix 2.0 4-shank only. You will need to train a new model with your own data if you have any other type of probe. (For example, Mouse2 from the figshare data is a Npix 1 dataset, you'll notice the trained model won't give good results on this mouse for that reason.)
