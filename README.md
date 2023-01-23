# UnitMatch
Ephys toolbox to match single units (electrophysiology) either within the same recording (oversplits) or across recordings

![image](https://user-images.githubusercontent.com/7621880/204564863-ab2e04f5-6d1c-4eb5-aa69-76eaa85cba0f.png)

This toolbox was created in November 2022-January 2023 by Enny H. van Beest with the purpose of creating a tool to match units across multiple recording(day)s, and/or merge oversplit units within the same day.
It uses spatial position and quality metrics similarity, both standard and unstandard metrics, to find possible matches. The alrogirthm is internally optimized using cross-correlation fingerprints (after the idea of Célian Bimbard).

This work is produced in the Carandini-Harris lab, in close collaboration with Célian Bimbard & Julie Fabre (with a close connection to her Bombcell repository).
Also Anna Lebedeva, Flora Takacs, Pip Coen, Reilly Tilbury Kenneth Harris and Matteo Carandini, as well as the rest of the laboratory are thanked for their feedback and contributions.
This work is supported by Marie Sklodowska-Curie 101022757 (Enny H. van Beest)

Toolboxes used for matching part:
- https://github.com/EnnyvanBeest/spikes (forked from https://github.com/cortex-lab/spikes, but also tested with the original spikes)
- https://github.com/kwikteam/npy-matlab
- https://github.com/Julie-Fabre/bombcell

Toolboxes used for other parts of pipeline:
- https://github.com/MouseLand/rastermap
- https://github.com/petersaj/AP_histology
- https://github.com/cortex-lab/allenCCF
- https://github.com/EnnyvanBeest/GeneralNeuropixelAnalysisFunctions

%% Usage
The UnitMatch function is called from within the function PrepareClusInfo, which will essentially read in cluster information, spike information (spikes toolbox), optionally runs Bombcells to perform quality metrics and find 'good' (single) units. 
You're welcome to use the pipeline as written by Enny (see MainAnalysis_Example), which will also produce some basic figures on unit distributions on the probe.
Alternatively, you can just input a list of Kilosort output directories to those recordings you want to find matches for (provided that the params file from kilosort output points to the correct raw data directory. If not please provide the RawDataPaths as well!)

