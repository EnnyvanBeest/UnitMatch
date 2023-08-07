# UnitMatch
Ephys toolbox to match single units (electrophysiology) either within the same recording (oversplits) or across recordings

![image](https://github.com/EnnyvanBeest/UnitMatch/blob/main/LogoAndExamples/Logo.svg)

This toolbox was created in November 2022-January 2023 by Enny H. van Beest as a tool to match units across multiple recording(day)s, and/or merge oversplit units within the same day. It uses spatial position and waveform-based parameters - both standard and unstandard metrics to the field - to find possible matches. The toolbox was further optimized between February-May 2023 in collaboration with Célian Bimbard.

We thank Julie Fabre, who implemented some changes to Bombcell - a toolbox for quality metrics and automated curation.
We also thank Anna Lebedeva, Flora Takacs, Pip Coen, Kenneth Harris, and Matteo Carandini, as well as the rest of the Carandini-Harris laboratory for their feedback and contributions.
This work was supported by Marie Sklodowska-Curie (101022757), Wellcome and NIH

Toolboxes used for matching part:
- https://github.com/EnnyvanBeest/spikes (forked from https://github.com/cortex-lab/spikes, but also tested with the original spikes)
- https://github.com/kwikteam/npy-matlab

Toolboxes that are very useful and integrated with UnitMatch:
- https://github.com/Julie-Fabre/bombcell

Toolboxes that are useful for other parts of analysis pipeline (but not necessary for just UnitMatch):
- https://github.com/MouseLand/rastermap
- https://github.com/petersaj/AP_histology
- https://github.com/cortex-lab/allenCCF
- https://github.com/EnnyvanBeest/GeneralNeuropixelAnalysisFunctions

### Usage

The UnitMatch algorithm is called with the function RunUnitMatch.m, with kilosort output directories and parameters as input. 
See the files `um_pipeline.m` (that calls`um_runUnitMatch.m`) or `RunUnitMatchAllDataPerMouse.m` for an example workflow. 


## Examples

Two recording sessions of same IMRO table. In the first recording this unit was oversplit, and UnitMatch will Merge it. Additonally it found it's match in the second recording.

![image](https://github.com/EnnyvanBeest/UnitMatch/blob/main/LogoAndExamples/Example1.bmp)

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

For commercial use please contact e.beest@ucl.ac.uk
