# GUI Reference Guide

## Shortcuts:
- [Right Arrow] - Move to NEXT unit pair in the slected list of matches (Unit A dropdown)
- [Left Arrow] - Move to the PREVIOUS unit pair in the slected list of matches (Unit A dropdown)
- [Up Arrow] - Move UP the Unit B probabilty list, e.g to a high probabilty match with Unit A
- [Down Arrow] - Move DOWN the Unit B probabilty list, e.g to a lower probabilty match with Unit A
- [q/m] - Set the selected pair of units as a MATCH
- [e/n] - Set the selected pair of units as a NON-MATCH
- [Enter/u] - Force update screen to show selected units 


## Unit Selection Guide:
To give maximum user control over the displayed infomaition the folowing can be used to select the unit pair.
- Unit A/B Session selection, you can select for each unit which session that unit is from using the dropwn down box.
- CV buttons, the are the 3 buttons which control what Cross-Verfiication pair the units are from the options are 'avg' '(1,2)' '(2,1)'
    avg, mean both Unit A and Unit B dispaly the average over their cv, (1,2) mean Unit A is in CV1 i.e the first half of that recording and
    Unit B is from CV2 i.e the second half of that recording. (2,1) is the reverse of (1,2)
- Unit A entry, this is a drop down box, the dropdown list is a list of all matches for the selected CV option and sessions
- Unit B entry, this is a drop down box , the dropdown list is a ordered list sorted from highest to lowest probabilty of being a matchwith Unit A

## Useful Infomation:
- The orginal Unit ID (as given in the list of good units) is displayed in the top right, the Unit ID in the unit selection is the UnitMatch ID
    this is so every unit has a unique ID
- Colors: Green = UnitA, Blue = UnitB
- Colors Histogram: ALL scores = Orange, Expected Matches = Magenta
- The raw waveforom plots are scaled sothe waveform plots have the same scale where the min and max values are those found in Unit A for the selected channels
- Stability, is the probabilty a unit matched with itself across CV. A value near 1 mean the unit mathces well with itself, while a value closer to 0 suggests the unit may not be stable.
- Hiding the raw waveform plots and UM score histograms, will speed up the refresh time of the GUI, slecting a toggle won't 
imediatley hide the plots, instead they will not update when toggled on.
