function LFP = LFPDepthFeatures(apFile, Opt)
% LFP = LFPDepthFeatures(apFile, Opt)
% LFP features per recording channel, for aligning histology to the probe:
%  - band power (delta, theta, gamma, ripple band; log10)
%  - ripple events per minute: narrowband 150-250 Hz events (envelope above median + Opt.RippleThrSD
%    robust SD for >= Opt.RippleMinMs, with at least Opt.RippleBandRatio times the power at 300-450 Hz,
%    so spike bursts (broadband) are not counted). Their rate peaks in CA1 stratum pyramidale.
%  - correlation with the nearest channel on the same shank (low at area boundaries)
% apFile: SpikeGLX .ap.bin/.ap.cbin. NP2.0 records wideband: the LFP is taken from the AP file.
%   NP1.0 (AP high-passed): the .lf.bin/.lf.cbin next to it is used (empty output if there is none).
% Opt.nSeg segments of Opt.SegSec s, spread over the recording, downsampled to 1 kHz.
% LFP: Shank (0-based), Depth (um, as the probe's site positions), Values (nChan x nFeatures), Names

if nargin < 2, Opt = struct; end
Def = struct('nSeg',10,'SegSec',20,'PieceSec',5,'RippleThrSD',5,'RippleMinMs',15,'RippleBandRatio',2);
for fn = fieldnames(Def)', if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end, end
LFP = [];
[fold,name,~] = fileparts(apFile);
meta = ReadMeta2(fold,'ap');
rawFile = apFile;
if isfield(meta,'imDatPrb_type') && ismember(str2double(meta.imDatPrb_type),[0 1020 1030 1100 1120 1121 1122 1123 1200 1300]) % NP1-like: separate LF file
    lf = [dir(fullfile(fold,[strrep(name,'.ap','.lf') '.cbin'])); dir(fullfile(fold,[strrep(name,'.ap','.lf') '.bin']))];
    if isempty(lf), warning('LFPDepthFeatures:NoLF','No LF file for %s',apFile); return; end
    rawFile = fullfile(lf(1).folder,lf(1).name);
end
% Channel positions (raw channel order, without the sync channel); shanks 250 um apart
chanPos = ChannelIMROConversion(fold,0,0);
nRec = size(chanPos,1);
Shank = floor(chanPos(:,1)/250); Depth = chanPos(:,2);
% nearest other channel on the same shank
nb = nan(nRec,1);
for c = 1:nRec
    same = find(Shank == Shank(c)); same(same == c) = [];
    if isempty(same), continue; end
    [~,k] = min(abs(Depth(same) - Depth(c)) + 0.01*abs(chanPos(same,1) - chanPos(c,1)));
    nb(c) = same(k);
end

[~,fs,nChan] = ReadRawChunk(rawFile,0,1);
if isfield(meta,'fileTimeSecs'), dur = str2double(meta.fileTimeSecs); else, dur = 600; end
fsOut = 1000; [p,q] = rat(fsOut/fs);
pieceLen = round(Opt.PieceSec*fs); nPiece = max(1,round(Opt.SegSec/Opt.PieceSec));
starts = round(linspace(0.05,0.95 - Opt.SegSec/max(dur,1),Opt.nSeg)*dur*fs);
[bR,aR] = butter(3,[150 250]/(fsOut/2),'bandpass');
[bH,aH] = butter(3,[300 450]/(fsOut/2),'bandpass');
minLen = round(Opt.RippleMinMs/1000*fsOut);
Pxx = 0; nP = 0; NbCorr = 0; nEvents = zeros(nRec,1); totSec = 0;
for s = 1:numel(starts)
    Y = cell(nPiece,1);
    for k = 1:nPiece
        X = ReadRawChunk(rawFile,starts(s) + (k-1)*pieceLen,pieceLen);
        X = double(X(:,1:min(nRec,nChan)));
        Y{k} = resample(X - mean(X,1),p,q); % anti-aliased downsampling to 1 kHz
    end
    Y = cat(1,Y{:});
    [P,f] = pwelch(Y,hann(fsOut),fsOut/2,fsOut,fsOut);
    Pxx = Pxx + P; nP = nP + 1;
    Yc = Y - median(Y,2); % common reference removed
    cc = nan(nRec,1);
    for c = find(~isnan(nb))', cc(c) = corr(Yc(:,c),Yc(:,nb(c))); end
    NbCorr = NbCorr + cc;
    % Ripple events
    envR = movmean(abs(hilbert(filtfilt(bR,aR,Yc))),round(0.008*fsOut));
    envH = movmean(abs(hilbert(filtfilt(bH,aH,Yc))),round(0.008*fsOut));
    thr = median(envR,1) + Opt.RippleThrSD*1.4826*mad(envR,1,1);
    for c = 1:nRec
        above = envR(:,c) > thr(c);
        d = diff([0; above; 0]); on = find(d == 1); off = find(d == -1) - 1;
        for e = find((off - on + 1) >= minLen)'
            [~,pk] = max(envR(on(e):off(e),c)); pk = pk + on(e) - 1;
            if envR(pk,c) >= Opt.RippleBandRatio*envH(pk,c), nEvents(c) = nEvents(c) + 1; end
        end
    end
    totSec = totSec + size(Y,1)/fsOut;
end
Pxx = Pxx/nP; NbCorr = NbCorr/numel(starts);
Bands = {'delta 1-4 Hz',[1 4]; 'theta 6-10 Hz',[6 10]; 'gamma 30-80 Hz',[30 80]; 'ripple band 150-250 Hz',[150 250]};
Values = nan(nRec,size(Bands,1)+2);
for b = 1:size(Bands,1)
    in = f >= Bands{b,2}(1) & f <= Bands{b,2}(2);
    Values(:,b) = log10(mean(Pxx(in,:),1))';
end
Values(:,end-1) = nEvents/totSec*60;
Values(:,end) = NbCorr;
LFP = struct('Shank',Shank,'Depth',Depth,'Values',Values,'Names',{[Bands(:,1)' {'ripple events per min','neighbour correlation'}]}, ...
    'File',rawFile,'Seconds',totSec);
end
