function Landmarks = TuningLandmarks(DataDir, Opt)
% Landmarks = TuningLandmarks(DataDir, Opt)
% Landmarks for the histology alignment (alignatlasdata_constrained, Opt.Landmarks) from the passive
% stimulation analyses of one probe on one day: RFData.mat (LoadSpikesRFBoot) and MPEPData.mat
% (Generalmpep_Analysis) in DataDir\<session>\ (DataDir = ...\Mouse\Date\ProbeX is also accepted:
% then the session folders of that date are searched for this probe).
% These analyses must not depend on the histology (they only use probe depths), so they can run
% before MainHistology.
% Each landmark: the fraction of units per depth with a property (Value 0/1 per unit, with its
% Shank (0-based) and Depth (um on the probe, as clusinfo.depth)), and the areas where that fraction
% should be high (AreaPaths: Allen structure id path tokens, e.g. '/669/' for visual cortex):
%   - RF: significant receptive field (PvalsAll < Opt.RFAlpha) -> visual areas
%   - visual / auditory responsiveness per mpep protocol: firing 0.05-0.5 s after stimulus onset
%     differs from the 0.2 s before (signed-rank over trials, p < Opt.RespAlpha)
%   - stimMergePureTonesAndFlickeringGratingsPlaids: locking to the grating temporal frequency
%     (tf_RayleighP) -> visual; to the sound modulation frequency (ModFreq_RayleighP) -> auditory
% Optogenetic protocols are skipped. Units are pooled over sessions per landmark type (a unit seen in
% several sessions counts once per session: fine for a depth profile).
if nargin < 2, Opt = struct; end
Def = struct('RFAlpha',0.05,'RespAlpha',0.01,'LockAlpha',1e-3,'MinUnits',20,'MinPositive',8,'Weight',0.2, ...
    'RespWin',[0.05 0.5],'BaseWin',[-0.2 0],'BinSec',0.01,'PreSec',0.2);
for fn = fieldnames(Def)'
    if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end
end
Visual   = {'/669/','/170/','/178/','/218/','/302/','/294/'}; % VIS, LGd, LGv, LP, SCs, SCm
Auditory = {'/247/','/475/','/4/'};                           % AUD, MG, IC

[dateDir,probe] = fileparts(DataDir);
if startsWith(probe,'Probe')
    rfFiles = dir(fullfile(dateDir,'*',probe,'RFData.mat'));
    mpFiles = dir(fullfile(dateDir,'*',probe,'MPEPData.mat'));
else
    rfFiles = dir(fullfile(DataDir,'**','RFData.mat'));
    mpFiles = dir(fullfile(DataDir,'**','MPEPData.mat'));
end

% Collected per landmark name: Shank, Depth, Value, AreaPaths, source files
Acc = struct('Name',{},'Shank',{},'Depth',{},'Value',{},'AreaPaths',{},'Files',{});
    function add(name,shank,depth,value,areas,file)
        ia = find(strcmp({Acc.Name},name),1);
        if isempty(ia)
            ia = numel(Acc) + 1;
            Acc(ia).Name = name; Acc(ia).AreaPaths = areas; Acc(ia).Shank = []; Acc(ia).Depth = []; Acc(ia).Value = []; Acc(ia).Files = {};
        end
        ok = ~isnan(value(:)) & ~isnan(depth(:));
        Acc(ia).Shank = [Acc(ia).Shank; double(shank(ok))];
        Acc(ia).Depth = [Acc(ia).Depth; double(depth(ok))];
        Acc(ia).Value = [Acc(ia).Value; double(value(ok))];
        Acc(ia).Files{end+1} = file;
    end

%% Receptive fields
for i = 1:numel(rfFiles)
    f = fullfile(rfFiles(i).folder,rfFiles(i).name);
    try
        r = load(f,'RFInfo'); R = r.RFInfo;
        p = getv(R,'PvalsAll'); if size(p,2) > 1, p = min(p,[],2); end
        add('RF (sig. receptive field)',col(getv(R,'Shank')),col(getv(R,'Depth')),double(p(:) < Opt.RFAlpha),Visual,f);
    catch ME
        warning('TuningLandmarks: could not use %s (%s)',f,ME.message)
    end
end

%% mpep protocols
for i = 1:numel(mpFiles)
    f = fullfile(mpFiles(i).folder,mpFiles(i).name);
    try
        m = load(f,'MpepInfo','Protocol');
        I = m.MpepInfo; x = m.Protocol.xfile; x = regexprep(x,'\.x$','');
        if contains(x,'Opti','IgnoreCase',true) || contains(x,'Laser','IgnoreCase',true), continue; end
        sh = col(getv(I,'Shank')); dp = col(getv(I,'depth'));
        if contains(x,'MergePureTonesAndFlickeringGratingsPlaids')
            % locking to the stimulus modulation: specific to the modality
            if hasv(I,'tf_RayleighP')
                add('mpep: locked to grating TF',sh,dp,double(min(getv(I,'tf_RayleighP'),[],2) < Opt.LockAlpha),Visual,f);
            end
            if hasv(I,'ModFreq_RayleighP')
                add('mpep: locked to sound modulation',sh,dp,double(min(getv(I,'ModFreq_RayleighP'),[],2) < Opt.LockAlpha),Auditory,f);
            end
            continue
        end
        if contains(x,'Grating') && ~contains(x,'Tones')
            name = 'mpep: grating responsive'; areas = Visual;
        elseif contains(x,'CosAmplModulatedSineWaves') || (contains(x,'Tone') && ~contains(x,'Grating'))
            name = 'mpep: sound responsive'; areas = Auditory;
        elseif contains(x,'Tones') && contains(x,'Grating')
            name = 'mpep: tone/grating responsive'; areas = [Visual Auditory]; % modality not separable here
        else
            continue
        end
        S = getv(I,'SpikesPerSec'); % units x trials x time
        if isempty(S) || ndims(S) ~= 3, continue; end
        t = ((1:size(S,3)) - 0.5)*Opt.BinSec - Opt.PreSec;
        base = mean(S(:,:,t >= Opt.BaseWin(1) & t < Opt.BaseWin(2)),3);
        resp = mean(S(:,:,t >= Opt.RespWin(1) & t < Opt.RespWin(2)),3);
        v = nan(size(S,1),1);
        for u = 1:size(S,1)
            ok = ~isnan(base(u,:)) & ~isnan(resp(u,:));
            if sum(ok) < 10, continue; end
            if all(resp(u,ok) == base(u,ok)), v(u) = 0; continue; end
            v(u) = double(signrank(resp(u,ok),base(u,ok)) < Opt.RespAlpha);
        end
        add(name,sh,dp,v,areas,f);
    catch ME
        warning('TuningLandmarks: could not use %s (%s)',f,ME.message)
    end
end

%% Output: landmarks with enough units
Landmarks = struct('Name',{},'Type',{},'Shank',{},'Depth',{},'Value',{},'AreaPaths',{},'AreaRegex',{},'Weight',{},'Files',{});
for k = 1:numel(Acc)
    if numel(Acc(k).Value) < Opt.MinUnits || sum(Acc(k).Value > 0) < Opt.MinPositive
        fprintf('Landmark %s not used: %d units, %d positive\n',Acc(k).Name,numel(Acc(k).Value),sum(Acc(k).Value > 0))
        continue
    end
    % Plain responsiveness (onset of a full-field stimulus) is found in many areas: half the weight of
    % the specific properties (receptive field, locking to the stimulus modulation)
    w = Opt.Weight; if endsWith(Acc(k).Name,'responsive'), w = w/2; end
    Landmarks(end+1) = struct('Name',Acc(k).Name,'Type','fraction','Shank',Acc(k).Shank,'Depth',Acc(k).Depth,'Value',Acc(k).Value, ...
        'AreaPaths',{Acc(k).AreaPaths},'AreaRegex','','Weight',w,'Files',{Acc(k).Files}); %#ok<AGROW>
    fprintf('Landmark %s: %d units, %d positive (%s)\n',Acc(k).Name,numel(Acc(k).Value),sum(Acc(k).Value),strjoin(unique(Acc(k).Files),', '))
end
end

function v = getv(S,name)
v = S.(name); % table variable or struct field
if iscell(v), v = cell2mat(v); end
v = double(v);
end

function tf = hasv(S,name)
if istable(S), tf = ismember(name,S.Properties.VariableNames); else, tf = isfield(S,name); end
end

function x = col(x)
x = x(:);
end
