function [Landmarks, CSD] = EvokedCSDLandmarks(DataDir, Opt)
% [Landmarks, CSD] = EvokedCSDLandmarks(DataDir, Opt)
% Visually evoked CSD of one probe on one day, as a landmark for the histology alignment.
% Per mpep session with a visual protocol (gratings, plaids): the CSD that Generalmpep_Analysis
% stored (MPEPData.mat, StimCSD), or computed here with StimulusCSD from the trial onsets on the
% probe clock (StimTiming). Per shank the session with the most trials is used.
% In visual cortex the earliest sink (just after the response onset) is in layer 4; in the superior
% colliculus in the superficial layers. Landmark: the early sink strength along the probe (sink
% positive, in baseline SDs), fitted to VIS layer 4 and the superficial SC (Allen acronyms,
% AreaRegex); only if the sink is clear (> Opt.ThrSD).
% Opt.Weight: 0.3 (within the cap on all landmarks); 0 shows it to the curator without fitting.
% DataDir: ...\Mouse\Date\ProbeX: MPEPData.mat of the sessions of that date are found in
%   ...\Mouse\Date\<session>\ProbeX\
% Landmarks: struct array (Name, Type 'profile', Shank, Depth, Value, AreaPaths, AreaRegex, Weight, Files)
% CSD: struct array per shank: Name, Shank, Depth, t, CSD, nTrials, OnsetMs, SinkDepth, SinkZ
if nargin < 2, Opt = struct; end
Def = struct('ThrSD',4,'Weight',0.3,'AreaRegex','^(VIS[a-z]*4|SCsg|SCzo|SCop)$');
for fn = fieldnames(Def)', if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end, end
Landmarks = struct('Name',{},'Type',{},'Shank',{},'Depth',{},'Value',{},'AreaPaths',{},'AreaRegex',{},'Weight',{},'Files',{});
CSD = struct('Name',{},'Shank',{},'Depth',{},'t',{},'CSD',{},'nTrials',{},'OnsetMs',{},'SinkDepth',{},'SinkZ',{});

[dateDir,probe] = fileparts(DataDir);
f = dir(fullfile(dateDir,'*',probe,'MPEPData.mat'));
All = []; Src = {};
for i = 1:numel(f)
    fn = fullfile(f(i).folder,f(i).name);
    w = who('-file',fn);
    if ~ismember('StimTiming',w), continue; end
    s = load(fn,'StimTiming'); S = s.StimTiming;
    if ~(contains(S.xfile,'Grating') || contains(S.xfile,'Plaid')), continue; end
    C = [];
    if ismember('StimCSD',w), s = load(fn,'StimCSD'); C = s.StimCSD; end
    if isempty(C)
        try
            C = StimulusCSD(S.RawFile,S.OnsetProbeSample);
        catch ME
            warning('EvokedCSDLandmarks:Failed','Visual CSD failed for %s: %s',fn,ME.message)
        end
    end
    for c = C(:)'
        if isempty(All), All = c; else, All(end+1) = c; end %#ok<AGROW>
        Src{end+1} = fn; %#ok<AGROW>
    end
end
if isempty(All), return; end
for sh = unique([All.Shank])
    k = find([All.Shank] == sh); [~,j] = max([All(k).nTrials]); c = All(k(j));
    CSD(end+1) = struct('Name',sprintf('visual CSD (%d trials)',c.nTrials),'Shank',sh,'Depth',c.Depth,'t',c.t, ...
        'CSD',c.CSD,'nTrials',c.nTrials,'OnsetMs',c.OnsetMs,'SinkDepth',c.SinkDepth,'SinkZ',c.SinkZ); %#ok<AGROW>
    fprintf('Visual CSD shank %d: %d trials, response onset %.0f ms, early sink %.1f SD at %.0f um (%s)\n', ...
        sh,c.nTrials,c.OnsetMs,c.SinkZ,c.SinkDepth,Src{k(j)})
    if c.SinkZ < Opt.ThrSD, continue; end
    ok = ~isnan(c.EarlySink);
    L = struct('Name','visual CSD early sink','Type','profile','Shank',repmat(sh,sum(ok),1),'Depth',c.Depth(ok), ...
        'Value',c.EarlySink(ok),'AreaPaths',{{}},'AreaRegex',Opt.AreaRegex,'Weight',Opt.Weight,'Files',{Src(k(j))});
    if isempty(Landmarks)
        Landmarks = L;
    else % one landmark over the shanks
        Landmarks.Shank = [Landmarks.Shank; L.Shank]; Landmarks.Depth = [Landmarks.Depth; L.Depth];
        Landmarks.Value = [Landmarks.Value; L.Value]; Landmarks.Files = [Landmarks.Files L.Files];
    end
end
end
