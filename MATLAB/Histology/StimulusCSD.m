function C = StimulusCSD(apFile, onsetSamples, Opt)
% C = StimulusCSD(apFile, onsetSamples, Opt)
% Stimulus-evoked LFP and CSD of one recording, per shank.
%  apFile: SpikeGLX .ap.bin/.cbin (NP1: the LF file next to it is read, see LFPSource)
%  onsetSamples: stimulus onsets in AP samples (e.g. StimTiming.OnsetProbeSample from Generalmpep_Analysis)
% The LFP (<Opt.LowPassHz, 1 kHz) is averaged around the onsets (at most Opt.MaxTrials, spread over
% the session; baseline Opt.BaseWin subtracted per trial) and the CSD computed (ComputeCSD).
% CSD per depth in its own baseline SDs (floor: the median), so ongoing activity that survives the
% average stays small. Response onset: first time (> 15 ms) the mean squared z over depth exceeds
% Opt.OnsetPowerX times its baseline level for 10 ms. Early sink: the largest sink in the first
% Opt.EarlyWinSec after that onset (Opt.SinkWin if no onset).
% C: struct array per shank: Shank, Depth, t, LFP, CSD (z), nTrials, OnsetMs, SinkWin, SinkDepth,
%   SinkZ, EarlySink (per depth, sink positive), File
if nargin < 3, Opt = struct; end
Def = struct('Win',[-0.1 0.3],'BaseWin',[-0.1 0],'SinkWin',[0.04 0.08],'OnsetPowerX',4,'EarlyWinSec',0.03, ...
    'LowPassHz',100,'MaxTrials',300,'MinTrials',20);
for fn = fieldnames(Def)', if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end, end
C = [];
if ~exist(apFile,'file') % .bin in the metadata, compressed .cbin on the server (or vice versa)
    alt = regexprep(apFile,'\.c?bin$',''); alt = [dir([alt '.cbin']); dir([alt '.bin'])];
    if isempty(alt), warning('StimulusCSD:NoRaw','No raw file %s',apFile); return; end
    apFile = fullfile(alt(1).folder,alt(1).name);
end
[rawFile, Shank, Depth, apToRaw] = LFPSource(apFile);
if isempty(rawFile), return; end
on = sort(onsetSamples(~isnan(onsetSamples)));
if numel(on) < Opt.MinTrials, return; end
if numel(on) > Opt.MaxTrials, on = on(round(linspace(1,numel(on),Opt.MaxTrials))); end
[~,fs,nChan] = ReadRawChunk(rawFile,0,1);
nRec = numel(Shank);
fsOut = 1000; [p,q] = rat(fsOut/round(fs)); % calibrated rates (e.g. 2500.01) give huge p/q; error negligible
[b,a] = butter(3,Opt.LowPassHz/(fsOut/2),'low');
pad = 0.2;
tOut = (round(Opt.Win(1)*fsOut):round(Opt.Win(2)*fsOut))/fsOut;
base = tOut >= Opt.BaseWin(1) & tOut < Opt.BaseWin(2);
Sum = zeros(nRec,numel(tOut)); n = 0;
for tr = 1:numel(on)
    s0 = round(on(tr)*apToRaw + (Opt.Win(1) - pad)*fs);
    nS = round((diff(Opt.Win) + 2*pad)*fs);
    if s0 < 0, continue; end
    X = ReadRawChunk(rawFile,s0,nS);
    if size(X,1) < nS, continue; end
    X = double(X(:,1:min(nRec,nChan)));
    Y = filtfilt(b,a,resample(X - mean(X,1),p,q));
    tY = (0:size(Y,1)-1)'/fsOut + Opt.Win(1) - pad;
    [~,i0] = min(abs(tY - tOut(1)));
    Y = Y(i0:i0+numel(tOut)-1,:);
    Y = Y - mean(Y(base,:),1);
    Sum = Sum + Y'; n = n + 1;
end
if n < Opt.MinTrials, return; end
CC = ComputeCSD(Sum/n,Shank,Depth,tOut,struct('BaseWin',Opt.BaseWin));
for c = 1:numel(CC)
    nz = std(CC(c).CSD(:,base),0,2,'omitnan'); nz = max(nz,median(nz,'omitnan'));
    Zc = CC(c).CSD./nz;
    r = mean(Zc.^2,1,'omitnan');
    above = r > Opt.OnsetPowerX*mean(r(base)) & tOut > 0.015;
    on1 = find(above & movmin(above,[0 9]),1);
    if isempty(on1), win = Opt.SinkWin; onMs = NaN; else, win = tOut(on1) + [0 Opt.EarlyWinSec]; onMs = 1000*tOut(on1); end
    sink = -mean(Zc(:,tOut >= win(1) & tOut < win(2)),2);
    [z,i] = max(sink);
    S = struct('Shank',CC(c).Shank,'Depth',CC(c).Depth,'t',tOut,'LFP',CC(c).LFP,'CSD',Zc,'nTrials',n,'OnsetMs',onMs, ...
        'SinkWin',win,'SinkDepth',CC(c).Depth(i),'SinkZ',z,'EarlySink',sink,'File',apFile);
    if isempty(C), C = S; else, C(end+1) = S; end %#ok<AGROW>
end
end
