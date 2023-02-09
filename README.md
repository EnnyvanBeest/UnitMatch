# UnitMatch
Ephys toolbox to match single units (electrophysiology) either within the same recording (oversplits) or across recordings

![image](https://github.com/EnnyvanBeest/UnitMatch/blob/main/Logo.png)

This toolbox was created in November 2022-January 2023 by Enny H. van Beest with the purpose of creating a tool to match units across multiple recording(day)s, and/or merge oversplit units within the same day. It uses spatial position and quality metrics similarity, both standard and unstandard metrics to the field, to find possible matches. The toolbox was further optimized in February 2023 in collaboration with Célian Bimbard - who sped up the code by a factor of almost infinity.

This work is produced in the Carandini-Harris lab, in close collaboration with Célian Bimbard & Julie Fabre (with a close connection to her Bombcell repository).
Also Anna Lebedeva, Flora Takacs, Pip Coen, Kenneth Harris, and Matteo Carandini, as well as the rest of the laboratory are thanked for their feedback and contributions.
This work is supported by Marie Sklodowska-Curie 101022757 (Enny H. van Beest)

Toolboxes used for matching part:
- https://github.com/EnnyvanBeest/spikes (forked from https://github.com/cortex-lab/spikes, but also tested with the original spikes)
- https://github.com/kwikteam/npy-matlab

Toolboxes very useful and integrated with UnitMatch:
- https://github.com/Julie-Fabre/bombcell

Toolboxes used for other parts of pipeline:
- https://github.com/MouseLand/rastermap
- https://github.com/petersaj/AP_histology
- https://github.com/cortex-lab/allenCCF
- https://github.com/EnnyvanBeest/GeneralNeuropixelAnalysisFunctions

%% Usage %%
The UnitMatch function is called from within the function PrepareClusInfo, which will essentially read in cluster information, spike information (spikes toolbox), optionally runs Bombcells to perform quality metrics and find 'good' (single) units. 
You're welcome to use the pipeline as written by Enny (see MainAnalysis_Example), which will also produce some basic figures on unit distributions on the probe.
Alternatively, you can just input a list of Kilosort output directories to those recordings you want to find matches for (provided that the params file from kilosort output points to the correct raw data directory. If not please provide the RawDataPaths as well!)

%% Examples %%
Two recording sessions of same IMRO table. In the first recording this unit was oversplit, and UnitMatch will Merge it. Additonally it found it's match in the second recording.

![image](https://github.com/EnnyvanBeest/UnitMatch/blob/main/Example1.png)
