# UnitMatch
Ephys toolbox to match single units (electrophysiology) either within the same recording (oversplits) or across recordings, by using spatial position and waveform-based parameters - both standard and unstandard metrics to the field. 

![image](https://github.com/EnnyvanBeest/UnitMatch/blob/main/MATLAB/LogoAndExamples/Logo.svg)

This toolbox was initially created between August 2022-January 2023 by Enny H. van Beest as a tool to match units across multiple recording(day)s, and/or merge oversplit units within the same day. The toolbox was further optimized between February-October 2023 in collaboration with Célian Bimbard.
A python version was created between November 2023 and January 2024 by Sam Dodgson. The python version only contains core UnitMatch code. Other functionality such as computing functional scores etc. for now is only available in the Matlab version.

We thank Julie Fabre, who implemented some changes to Bombcell - a toolbox for quality metrics and automated curation.
We also thank Anna Lebedeva, Flora Takacs, Pip Coen, Kenneth Harris, and Matteo Carandini, as well as the rest of the Carandini-Harris laboratory for their feedback and contributions.

For this work, van Beest was supported by a Marie Sklodowska-Curie fellowship (van Beest, 101022757), Bimbard by a European Molecular Biology Organization (ALTF 740-2019), and Carandini and Harris by a Wellcome Investigator Award.


Preprint: https://www.biorxiv.org/content/10.1101/2023.10.12.562040v1

Video: https://www.youtube.com/watch?v=4c_dgTZcBaQ&list=PLfhWmWntvjl7kljKozClpjS29DoY8V5pB&index=23&t=21s

Below instructions are for the Matlab version. Please see the python folder for more information on the Python version.
### Dependencies on other toolboxes/repositories
Toolboxes used for matching part:
- https://github.com/kwikteam/npy-matlab

Toolboxes that are very useful and integrated with UnitMatch:
- https://github.com/Julie-Fabre/bombcell
- https://github.com/EnnyvanBeest/spikes (forked from https://github.com/cortex-lab/spikes, but also tested with the original spikes)

Toolboxes that are useful for other parts of analysis pipeline (but not necessary for just UnitMatch):
- https://github.com/MouseLand/rastermap
- https://github.com/petersaj/AP_histology
- https://github.com/cortex-lab/allenCCF
- https://github.com/EnnyvanBeest/GeneralNeuropixelAnalysisFunctions

### Usage
We included a DEMO_UNITMATCH.m script, which hopefully clarifies how to use UnitMatch and smoothly integrate it in your existing analysis pipeline. UnitMatch requires at minimum a cell array with path names, pointing do the different recordings for which you want to try and track units. Each of these paths (typically a Kilosort output folder) should contain a subfolder called 'RawWaveforms' in which for every unit (/cluster) there is a NPY file containing two average waveforms per recording site (spikewidth X nRecordingSites X 2). The two average waveforms are preferably from the first versus second half of a recording. Additionally UnitMatch needs a 'clusinfo' struct with the fields cluster_id, Good_ID (whether to include it or not), RecSesID (which recording session), and Probe (can be a vector of ones if just one probe was used). It also needs some parameters, which can be filled in using DefaultParametersUnitMatch.m. 

Also using spikeGLX and Neuropixels, and spikesorting your data with (some form of) Kilosort? You are lucky!
We included two example pipelines (see the folder 'ExampleAnalysisPipelines'), to show how we go from raw data (SpikeGLX + Neuropixels) to Kilosorted data ((Py)Kilosort) to curated single units (Bombcell), from where we can smoothly run UnitMatch.

### Output
UnitMatch has two main outputs:
1. A matching table. This contains for every included pair of units each of the similarity scores, the probability of being a match, and some extra information (e.g. functional scores, if using our functional score validation analysis).
2. A UniqueIDConversion, which is a struct containing the original cluster identities (as was defined in clusinfo.cluster_id), the 'UniqueID' which - if using the 'AssignUniqueID' function - gives matching clusters across recording days the same UniqueID. Useful for later analysis! It contains extra information that might be useful, such as whether it was considered a 'good unit', and what recording session it was found in.

### Modules
After the initial UnitMatch algorithm, you can use different 'Modules' to check and validate UnitMatch' output. Each module can be used simply by entering the path to the UnitMatch output. N.B. for some of the modules you need functional data in the format of (Py)Kilosort output or need to have used Bombcell.
1. AssignUniqueID (see Output.2)
2. EvaluatingUnitMatch: Within session cross-validation, or comparison to using Kilosort in the stitched/concatenated way
3. QualityMetricsROCs: Check for every quality metric from Bombcell whether this affected a unit's chance of being a match with another unit
4. ComputeFunctionalScores: Computes functional scores (autocorrelograms, reference population correlations, firing rate differences, if available natural images responses) similarities for each pair of units. Will be added to the matching table output.
5. DrawPairsUnitMatch: Draws and saves a figure for every matching pair that was found. It will also include some 'doubt cases', such as pairs that had very high functional scores, yet were not found to be a match.
6. FigureFlick: Useful for manual curation. A column will be added to the Matching table with the user's judgment for pairs of units. Only works after running the DrawPairsUnitMatch module.

### Examples

Two recording sessions of same IMRO table. In the first recording this unit was oversplit, and UnitMatch will Merge it. Additonally it found it's match in the second recording.

![image](https://github.com/EnnyvanBeest/UnitMatch/blob/main/MATLAB/LogoAndExamples/Example1.bmp)

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

For commercial use please contact e.beest@ucl.ac.uk

