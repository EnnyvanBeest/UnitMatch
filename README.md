# UnitMatch
Ephys toolbox to match single units (electrophysiology) either within the same day (over splits) or across two days

![image](https://user-images.githubusercontent.com/7621880/204564863-ab2e04f5-6d1c-4eb5-aa69-76eaa85cba0f.png)

This toolbox was created in November 2022 by Enny H. van Beest, with the purpose of creating a tool to match units across multiple days, and/or merge oversplit units within the same day.
It uses spatial position and quality metrics similarity to find possible matches. The alrogirthm is internally optimized using cross-correlation fingerprints.

This work is produced in the Carandini-Harris lab, in close collaboration with Célian Bimbard & Julie Fabre.
Also Anna Lebedeva, Flora Takacs, Pip Coen, Reilly Tilbury Kenneth Harris and Matteo Carandini, as well as the rest of the laboratory are thanked for their feedback and contributions.
This work is supported by Marie Sklodowska-Curie 101022757 (Enny H. van Beest)

Toolboxes used for matching part:
- https://github.com/EnnyvanBeest/spikes (forked from https://github.com/cortex-lab/spikes)
- https://github.com/kwikteam/npy-matlab
- https://github.com/Julie-Fabre/bombcell

Toolboxes used for other parts of pipeline:
- https://github.com/MouseLand/rastermap
- https://github.com/petersaj/AP_histology
- https://github.com/cortex-lab/allenCCF
- https://github.com/EnnyvanBeest/GeneralNeuropixelAnalysisFunctions
