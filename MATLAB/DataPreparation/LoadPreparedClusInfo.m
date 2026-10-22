function [clusinfo, sp, Params, recsesidxIncluded] = LoadPreparedClusInfo(KiloSortPaths,Params)

%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1,length(KiloSortPaths));
addthis=0;
StoreRecSesID = cell(1,length(KiloSortPaths));
for subsesid=1:length(KiloSortPaths)
     if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy'))) || ~exist(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'))
        continue
    end

    disp(['Loading clusinfo for ' KiloSortPaths{subsesid}])
    tmp = load(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'clusinfo');
    clusinfo{subsesid} = tmp.clusinfo;

    % Replace recsesid with subsesid
    clusinfo{subsesid}.RecSesID = clusinfo{subsesid}.RecSesID+addthis;
    StoreRecSesID{subsesid} = unique(clusinfo{subsesid}.RecSesID);
    addthis=max(clusinfo{subsesid}.RecSesID);
end



% Add all cluster information in one 'cluster' struct - can be used for further analysis
clusinfo = [clusinfo{:}];
if isempty(clusinfo)
    clusinfo = [];
    sp = [];
    return
end
clusinfoNew = struct;
fields = fieldnames(clusinfo(1));
for fieldid=1:length(fields)
    try
        eval(['clusinfoNew.' fields{fieldid} '= cat(1,clusinfo(:).' fields{fieldid} ');'])
    catch ME
        try
            eval(['clusinfoNew.' fields{fieldid} '= cat(2,clusinfo(:).' fields{fieldid} ');'])
        catch ME
            keyboard
        end
    end
end
clusinfo = clusinfoNew;
clear clusinfoNew

%% Here we're going to actually load in all the sessions requested - sp
sp = cell(1,length(KiloSortPaths));
countid=1;
for subsesid=1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))  || ~exist(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'))
        continue
    end
    disp(['Loading spike data for ' KiloSortPaths{subsesid}])
    tmp = load(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'sp');
    sp{subsesid} = tmp.sp;
    % Replace recsesid with subsesid
    sp{subsesid}.RecSes = repmat(countid,size(sp{subsesid}.RecSes));
    countid=countid+1;
end
% Add all spikedata in one spikes struct - can be used for further analysis
sp = [sp{:}];
spnew = struct;

fields = fieldnames(sp(1));
for fieldid=1:length(fields)
    try
        eval(['spnew.' fields{fieldid} '= cat(1,sp(:).' fields{fieldid} ');'])
    catch ME
        if strcmp(ME.message,'Out of memory.')
            eval(['spnew.' fields{fieldid} ' = sp(1).' fields{fieldid} ';'])
            for tmpid = 2:length(sp)
                eval(['spnew.' fields{fieldid} ' = cat(1,spnew.' fields{fieldid} ', sp(tmpid).' fields{fieldid} ');'])
            end
        else
            eval(['spnew.' fields{fieldid} '= cat(2,sp(:).' fields{fieldid} ');'])
        end
    end
end
sp = spnew;
sp.sample_rate = sp.sample_rate(1);
clear spnew

% Save out sp unique cluster
recsesidxIncluded = false(1,length(KiloSortPaths));
nclus = length(clusinfo.cluster_id);
RawPathsUsed = {};
sp.UniqClu = sp.clu;
sp.SpikeOnProbe = nan(size(sp.clu));
clusinfo.UniqueID = (1:length(clusinfo.cluster_id))';
clusinfo.IMROID = repmat(0,length(clusinfo.cluster_id),1);
clusinfo.AdjustedDepth = clusinfo.depth;
MaxDriftUm = 300; % drift estimates larger than this (um) are not trusted: depth is then left uncorrected
SN = []; % default; overwritten inside UnitMatch block below
if Params.UnitMatch
    disp('Assigning correct unique ID')
    PartsPath = strsplit(Params.SaveDir,'\');
    UMOutputAll = dir(fullfile(PartsPath{1},PartsPath{2},PartsPath{3},'**','UnitMatch','UnitMatch.mat'));
    nAssigned = 0;
    for imroid = 1:length(UMOutputAll)
        if nAssigned >= numel(KiloSortPaths) % Already have all information
            break
        end
        IMROId = strsplit(UMOutputAll(imroid).folder,'\');
        IMROId = strsplit(IMROId{end-1},'_');
        if strcmp(IMROId,'AllIMRO') & ~Params.separateIMRO
            IMROId = 0;
        else
            IMROId = str2num(IMROId{end});
        end
        if Params.separateIMRO & isempty(IMROId)
            continue
        end

        UMOutput = load(fullfile(UMOutputAll(imroid).folder,UMOutputAll(imroid).name),'UMparam');
        UMparam = UMOutput.UMparam;
        SN = UMparam.AllProbeSN;
        recsesidx = find(ismember(UMparam.KSDir,KiloSortPaths)); % Find which recses id they should have in UM output
        recsesidx2 = find(ismember(KiloSortPaths,UMparam.KSDir));
        if ~any(recsesidx)
            continue
        else
            nAssigned = nAssigned + 1;
        end

        UMOutput = load(fullfile(UMOutputAll(imroid).folder,UMOutputAll(imroid).name),'UMparam','UniqueIDConversion');       
        UniqueIDConversion = UMOutput.UniqueIDConversion;
        
        TheseClus = find(ismember(clusinfo.RecSesID,[StoreRecSesID{recsesidx2}]));
        if ~isempty(TheseClus)
            clusinfo.UniqueID(TheseClus) = UniqueIDConversion.UniqueID(ismember(UniqueIDConversion.recsesAll,recsesidx))'; %Assign correct UniqueID
            clusinfo.IMROID(TheseClus) = repmat(IMROId,sum(ismember(clusinfo.RecSesID,[StoreRecSesID{recsesidx2}])),1);
            % Depth relative to the first recording of this UnitMatch run (Z = depth axis only)
            DriftZ = DriftPerRecording(UMparam,MaxDriftUm);
            for rrid = 1:numel(recsesidx)
                TheseClus = find(ismember(clusinfo.RecSesID,[StoreRecSesID{recsesidx2(rrid)}]));
                clusinfo.AdjustedDepth(TheseClus) = clusinfo.depth(TheseClus) + DriftZ(recsesidx(rrid));
            end
        end
        sp.AdjustedspikeDepth = sp.spikeDepths;
        for clusid=1:length(TheseClus)
            sp.UniqClu(sp.clu==clusinfo.cluster_id(TheseClus(clusid)) & sp.RecSes==clusinfo.RecSesID(TheseClus(clusid))) = clusinfo.UniqueID(TheseClus(clusid));
            sp.SpikeOnProbe(sp.clu==clusinfo.cluster_id(TheseClus(clusid)) & sp.RecSes==clusinfo.RecSesID(TheseClus(clusid))) = clusinfo.ProbeID(TheseClus(clusid));
            sp.AdjustedspikeDepth(sp.clu==clusinfo.cluster_id(TheseClus(clusid)) & sp.RecSes==clusinfo.RecSesID(TheseClus(clusid))) = clusinfo.AdjustedDepth(TheseClus(clusid));

        end
        recsesidxIncluded(recsesidx2) = 1;
        for rrid = 1:numel(recsesidx)
            if isstruct(UMparam.RawDataPaths{recsesidx(rrid)})
                UMparam.RawDataPaths{recsesidx(rrid)} = fullfile(UMparam.RawDataPaths{recsesidx(rrid)}.folder,UMparam.RawDataPaths{recsesidx(rrid)}.name);
            end
        end

        RawPathsUsed = {RawPathsUsed{:} UMparam.RawDataPaths{recsesidx}};
    end
end
Params.RawDataPaths = RawPathsUsed;
Params.KSDir = KiloSortPaths;
Params.SN = SN;
%% add nans for missing UID values to not confuse them with other UIDs
if any(recsesidxIncluded==0)
    clusinfo.UniqueID(ismember(clusinfo.RecSesID,find(recsesidxIncluded==0))) = nan;
    clusinfo.IMROID(ismember(clusinfo.RecSesID,find(recsesidxIncluded==0))) = nan;
else % Organize units; only one entry per UID to avoid double counts of neurons
    [UOpt,ID1,ID2] = unique(cat(2,clusinfo.UniqueID,clusinfo.IMROID),'rows','stable'); % Look for unique ID/IMRO pairs

    thesefields = fieldnames(clusinfo);
    for fid = 1:numel(thesefields)
        clusinfo.(thesefields{fid}) =  clusinfo.(thesefields{fid})(ID1);
    end
end




return

function DriftZ = DriftPerRecording(UMparam,MaxDriftUm)
% Depth (Z) offset of every recording of a UnitMatch run relative to its first recording, in um
% (add to a unit's depth). UMparam.drift is (recording pair k -> k+1, [X Y Z], pass):
%  pass 1: correction applied to recording k+1, estimated against recording k after k's own
%          correction, i.e. relative to the first recording (if the chain is not broken)
%  pass 2: what remains between recordings k and k+1 after pass 1 (measured, not applied)
% Offset of recording r: its pass-1 correction. Where pass 1 had too few matches for a pair (NaN:
% recording k+1 not corrected, and later pass-1 corrections are relative to it), the pass-2 estimate
% of that pair bridges the gap and is added to all later recordings. Pass-2 remainders of other
% pairs are not added (they are mostly noise and would accumulate over the chain).
% Pairs without any estimate count as 0 (a warning is given).
nRec = numel(UMparam.KSDir);
DriftZ = zeros(nRec,1);
if ~isfield(UMparam,'drift') || isempty(UMparam.drift) || nRec < 2
    return
end
d = UMparam.drift;
nPair = min(nRec-1,size(d,1));
Applied = d(1:nPair,3,1);
if size(d,3) >= 2, Rest = d(1:nPair,3,end); else, Rest = zeros(nPair,1); end
Broken = isnan(Applied) & isnan(Rest);
if any(Broken)
    warning('LoadPreparedClusInfo:DriftChainBroken', ...
        'No drift estimate between recordings %s and the next (too few matches): later recordings may not be aligned to the first', ...
        mat2str(find(Broken)'))
end
Bridge = Rest .* isnan(Applied); % pass-2 estimate only for pairs that pass 1 could not align
Bridge(isnan(Bridge)) = 0; Applied(isnan(Applied)) = 0;
DriftZ(2:nPair+1) = Applied + cumsum(Bridge);
TooLarge = abs(DriftZ) > MaxDriftUm;
if any(TooLarge)
    warning('LoadPreparedClusInfo:DriftTooLarge', ...
        'Drift estimate of recording(s) %s larger than %d um (%s um): not used, depth left uncorrected', ...
        mat2str(find(TooLarge)'),MaxDriftUm,mat2str(round(DriftZ(TooLarge)')))
    DriftZ(TooLarge) = 0;
end