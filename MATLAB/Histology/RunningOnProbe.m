function Run = RunningOnProbe(DataDir, apFile)
% Run = RunningOnProbe(DataDir, apFile)
% Running speed on the probe clock for one recording, from all its sessions (VR, RF mapping, mpep,
% ...; each with its own Timeline): RunTiming.mat in DataDir's date folder, saved by SaveRunTiming
% (MainAnalysis_SpatialNavigation, run after the sync in MatchSyncedRecording).
% DataDir: ...\Mouse\Date\ProbeX (sessions in ...\Mouse\Date\<session>\ProbeX\); apFile: the raw
% AP file (sessions recorded into another file are not used).
% Run: struct with ProbeSample (AP samples) and Speed (cm/s), sorted; empty if none
Run = [];
[dateDir,probe] = fileparts(DataDir);
f = dir(fullfile(dateDir,'*',probe,'RunTiming.mat'));
[~,want] = fileparts(regexprep(apFile,'\.c?bin$',''));
ps = []; v = [];
for i = 1:numel(f)
    s = load(fullfile(f(i).folder,f(i).name),'RunTiming'); R = s.RunTiming;
    [~,have] = fileparts(regexprep(R.RawFile,'\.c?bin$',''));
    if ~strcmp(have,want) || isempty(R.Speed), continue; end
    ps = [ps; R.ProbeSample(:)]; v = [v; R.Speed(:)]; %#ok<AGROW>
end
if isempty(ps), return; end
[ps,o] = sort(ps);
Run = struct('ProbeSample',ps,'Speed',v(o));
end
