function LFP = LFPDepthFeatures(apFile, Opt)
% LFP = LFPDepthFeatures(apFile, Opt)
% LFP features per recording channel, for aligning histology to the probe:
%  - band power (delta, theta, gamma, ripple band; log10)
%  - ripple events per minute: narrowband 150-250 Hz events (envelope above median + Opt.RippleThrSD
%    robust SD for >= Opt.RippleMinMs, with at least Opt.RippleBandRatio times the power at 300-450 Hz,
%    so spike bursts (broadband) are not counted). Their rate peaks in CA1 stratum pyramidale.
%  - correlation with the nearest channel on the same shank (low at area boundaries)
%  - ripple-triggered CSD per shank (LFP.RippleCSD): the sharp wave (1-40 Hz) averaged around the
%    ripple events of the channel with the most events. In CA1 the sharp wave has a sink in stratum
%    radiatum, below a source at the pyramidal layer (where the ripples are): Valid if both are
%    there (> Opt.CSDThrSD x noise). Ripple-band events without this pattern (e.g. thalamic bursts)
%    are not CA1.
%  - theta-trough-triggered CSD per shank (LFP.ThetaCSD)
%  - dentate spikes per minute and the dentate-spike CSD per shank (LFP.DentateCSD): DG hilus
%  - slow-wave (delta trough) CSD per shank (LFP.SlowCSD): deep-layer sink in the neocortex
% Brain state: ripples (and sharp waves) occur mostly during quiet rest, theta during running. The
% Opt.nSeg short segments are ranked by the theta/delta power ratio on the channel with the most theta
% (a reliable correlate of locomotion; no behavioural sync needed): the theta CSD uses the most
% active Opt.StateFrac of the segments, ripples and the ripple CSD all others; band powers use all.
% With the running speed on the probe clock (Opt.Run: ProbeSample in AP samples, Speed in cm/s; see
% RunningOnProbe) a third of the segments is taken from running (mean speed > Opt.RunThr) and a third
% from still periods (max speed < Opt.StillThr), the rest spread over the recording; then running
% segments are the active state and still segments the quiet state.
% apFile: SpikeGLX .ap.bin/.ap.cbin. NP2.0 records wideband: the LFP is taken from the AP file.
%   NP1.0 (AP high-passed): the .lf.bin/.lf.cbin next to it is used (empty output if there is none).
% Opt.nSeg segments of Opt.SegSec s, spread over the recording, downsampled to 1 kHz.
% LFP: Shank (0-based), Depth (um, as the probe's site positions), Values (nChan x nFeatures), Names,
%   RippleCSD (struct per shank: ComputeCSD output + RefDepth, nEvents, SinkDepth, SourceDepth,
%   SinkZ, SourceZ, Valid), ThetaCSD, State (theta/delta ratio per segment, Quiet, Active, RefDepth)

if nargin < 2, Opt = struct; end
Def = struct('nSeg',90,'SegSec',4,'PieceSec',4,'StateFrac',1/3,'Run',[],'RunThr',5,'StillThr',1,'RippleThrSD',5,'RippleMinMs',15,'RippleBandRatio',2, ...
    'RippleNbUm',40,'RippleNbSD',3,'RippleMinNbFrac',0.5,'CSDMinEvents',25,'CSDThrSD',3,'CSDWinSec',0.15, ...
    'DSThrSD',6,'DSWidthMs',[5 25],'DSMinEvents',20,'DSCSDThrSD',5,'SlowExclThetaSD',10,'SlowHippoThetaSD',20,'SlowCSDThrSD',5);
for fn = fieldnames(Def)', if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end, end
LFP = [];
[fold,~,~] = fileparts(apFile);
meta = ReadMeta2(fold,'ap');
[rawFile, Shank, Depth, apToRaw] = LFPSource(apFile);
if isempty(rawFile), return; end
chanPos = ChannelIMROConversion(fold,0,0);
nRec = numel(Shank);
% nearest other channel on the same shank
nb = nan(nRec,1);
for c = 1:nRec
    same = find(Shank == Shank(c)); same(same == c) = [];
    if isempty(same), continue; end
    [~,k] = min(abs(Depth(same) - Depth(c)) + 0.01*abs(chanPos(same,1) - chanPos(c,1)));
    nb(c) = same(k);
end
% channels within Opt.RippleNbUm on the same shank (for the ripple consensus)
nbR = arrayfun(@(c) setdiff(find(Shank == Shank(c) & abs(Depth - Depth(c)) <= Opt.RippleNbUm),c),(1:nRec)','Uni',0);

[~,fs,nChan] = ReadRawChunk(rawFile,0,1);
if isfield(meta,'fileTimeSecs'), dur = str2double(meta.fileTimeSecs); else, dur = 600; end
fsOut = 1000; [p,q] = rat(fsOut/round(fs)); % calibrated rates (e.g. 2500.01) give huge p/q; error negligible
pieceLen = round(Opt.PieceSec*fs); nPiece = max(1,round(Opt.SegSec/Opt.PieceSec));
starts = round(linspace(0.05,0.95 - Opt.SegSec/max(dur,1),Opt.nSeg)*dur*fs);
segLabel = zeros(numel(starts),1); % behaviour: 1 running, -1 still, 0 not selected by behaviour
if ~isempty(Opt.Run) && numel(Opt.Run.Speed) > 100
    rs = double(Opt.Run.ProbeSample(:))*apToRaw; v = double(Opt.Run.Speed(:));
    segLen = round(Opt.SegSec*fs);
    cand = (rs(1):round(segLen/2):rs(end) - segLen)'; % half-overlapping candidate windows
    mv = arrayfun(@(c) mean(v(rs >= c & rs < c + segLen)),cand);
    xv = arrayfun(@(c) max([v(rs >= c & rs < c + segLen); Inf*isempty(v(rs >= c & rs < c + segLen))]),cand);
    runC = cand(mv > Opt.RunThr); stillC = cand(xv < Opt.StillThr);
    nPer = round(Opt.nSeg*Opt.StateFrac);
    if numel(runC) >= 10 && numel(stillC) >= 10
        pick = @(c) c(unique(round(linspace(1,numel(c),min(nPer,numel(c))))));
        sR = pick(runC); sS = pick(stillC); nU = Opt.nSeg - numel(sR) - numel(sS);
        sU = round(linspace(0.05,0.95 - Opt.SegSec/max(dur,1),nU)*dur*fs)';
        starts = round([sR; sS; sU])'; segLabel = [ones(numel(sR),1); -ones(numel(sS),1); zeros(nU,1)];
    end
end
[bR,aR] = butter(3,[150 250]/(fsOut/2),'bandpass');
[bH,aH] = butter(3,[300 450]/(fsOut/2),'bandpass');
minLen = round(Opt.RippleMinMs/1000*fsOut);
[bS,aS] = butter(3,[1 40]/(fsOut/2),'bandpass'); dsS = 4; % sharp wave, kept at 250 Hz for the CSD
Pxx = 0; nP = 0; NbCorr = 0; totSec = 0;
Ylow = cell(numel(starts),1); evPk = cell(numel(starts),1);
segTheta = zeros(numel(starts),nRec); segDelta = zeros(numel(starts),nRec); segSec = zeros(numel(starts),1);
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
    segTheta(s,:) = mean(P(f >= 6 & f <= 10,:),1); segDelta(s,:) = mean(P(f >= 1 & f <= 4,:),1);
    Yc = Y - median(Y,2); % common reference removed
    cc = nan(nRec,1);
    for c = find(~isnan(nb))', cc(c) = corr(Yc(:,c),Yc(:,nb(c))); end
    NbCorr = NbCorr + cc;
    % Ripple events
    envR = movmean(abs(hilbert(filtfilt(bR,aR,Yc))),round(0.008*fsOut));
    envH = movmean(abs(hilbert(filtfilt(bH,aH,Yc))),round(0.008*fsOut));
    thr = median(envR,1) + Opt.RippleThrSD*1.4826*mad(envR,1,1);
    thrNb = median(envR,1) + Opt.RippleNbSD*1.4826*mad(envR,1,1);
    evPk{s} = cell(nRec,1);
    for c = 1:nRec
        above = envR(:,c) > thr(c);
        d = diff([0; above; 0]); on = find(d == 1); off = find(d == -1) - 1;
        for e = find((off - on + 1) >= minLen)'
            [~,pk] = max(envR(on(e):off(e),c)); pk = pk + on(e) - 1;
            % a ripple spans the pyramidal layer: also on >= Opt.RippleMinNbFrac of the channels within Opt.RippleNbUm
            % (channel-specific high-frequency bursts are not ripples)
            if envR(pk,c) >= Opt.RippleBandRatio*envH(pk,c) && ~isempty(nbR{c}) && mean(envR(pk,nbR{c}) > thrNb(nbR{c})) >= Opt.RippleMinNbFrac
                evPk{s}{c}(end+1) = pk;
            end
        end
    end
    Ylow{s} = single(filtfilt(bS,aS,Y)); Ylow{s} = Ylow{s}(1:dsS:end,:);
    segSec(s) = size(Y,1)/fsOut; totSec = totSec + segSec(s);
end

%% Brain state per segment: theta/delta ratio on the channel with the most theta power
[~,refT] = max(sum(segTheta,1));
ratio = segTheta(:,refT)./max(segDelta(:,refT),eps);
if any(segLabel == 1) && any(segLabel == -1) % running speed known
    Active = segLabel == 1; Quiet = segLabel == -1; src = 'running speed';
else
    Active = ratio >= prctile(ratio,100*(1 - Opt.StateFrac));
    Quiet = ~Active; % ripples: all but the theta (running) segments, to keep enough events
    src = 'theta/delta';
end
State = struct('Source',src,'ThetaDeltaRatio',ratio,'Quiet',Quiet,'Active',Active,'RefDepth',Depth(refT));
% ripple events in the quiet segments
nEvents = zeros(nRec,1);
for s = find(Quiet)', nEvents = nEvents + cellfun(@numel,evPk{s}); end
quietSec = sum(segSec(Quiet));

%% Ripple-triggered CSD per shank (events of the channel with the most events)
RippleCSD = [];
hw = round(Opt.CSDWinSec*fsOut/dsS); tS = (-hw:hw)/(fsOut/dsS);
for sh = unique(Shank(:))'
    ch = find(Shank == sh);
    [nRef,k] = max(nEvents(ch)); ref = ch(k);
    if nRef < Opt.CSDMinEvents, continue; end
    Sum = zeros(numel(ch),numel(tS)); n = 0;
    for s = find(Quiet)'
        for pk = round(evPk{s}{ref}/dsS)
            if pk - hw < 1 || pk + hw > size(Ylow{s},1), continue; end
            Sum = Sum + double(Ylow{s}(pk-hw:pk+hw,ch))'; n = n + 1;
        end
    end
    if n < Opt.CSDMinEvents, continue; end
    C = ComputeCSD(Sum/n,Shank(ch),Depth(ch),tS,struct('BaseWin',[-Inf -0.08]));
    if isempty(C), continue; end
    % Sharp wave at the ripple (+-20 ms): source at the pyramidal layer (ripple channel), sink in
    % stratum radiatum 50-450 um deeper (lower depth on the probe)
    m = mean(C.CSD(:,abs(tS) <= 0.02),2)/max(C.Noise,eps);
    rd = Depth(ref);
    inSrc = C.Depth >= rd - 100 & C.Depth <= rd + 150;
    inSnk = C.Depth >= rd - 450 & C.Depth <= rd - 50;
    C.RefDepth = rd; C.nEvents = n;
    C.SourceZ = NaN; C.SourceDepth = NaN; C.SinkZ = NaN; C.SinkDepth = NaN;
    if any(inSrc & ~isnan(m)), [C.SourceZ,i] = max(m .* (inSrc./inSrc)); C.SourceDepth = C.Depth(i); end
    if any(inSnk & ~isnan(m)), [C.SinkZ,i] = min(m .* (inSnk./inSnk)); C.SinkDepth = C.Depth(i); end
    C.Valid = C.SourceZ >= Opt.CSDThrSD && C.SinkZ <= -Opt.CSDThrSD;
    if isempty(RippleCSD), RippleCSD = C; else, RippleCSD(end+1) = C; end %#ok<AGROW>
end

%% Dentate spikes (quiet segments): brief (Opt.DSWidthMs), large (> Opt.DSThrSD robust SD) positive
% LFP events, largest in the DG hilus, with sinks in the molecular layer; not within 50 ms of a
% ripple (the positive sharp-wave deflection at the CA1 pyramidal layer is not a dentate spike).
% Per shank on the channel with the most of them; Valid if the dentate-spike CSD has a source at that
% channel (+-60 um) and a sink 60-350 um above it (molecular layer of the upper blade).
DentateCSD = []; nDS = zeros(nRec,1);
fsL = fsOut/dsS;
dsPk = cell(numel(starts),nRec);
for s = find(Quiet)'
    nearRip = false(size(Ylow{s},1),1);
    for c = 1:nRec % +-50 ms around any ripple
        for pk = round(evPk{s}{c}/dsS), nearRip(max(1,pk-round(0.05*fsL)):min(end,pk+round(0.05*fsL))) = true; end
    end
    Zs = (double(Ylow{s}) - median(Ylow{s},1))./(1.4826*mad(double(Ylow{s}),1,1));
    for c = 1:nRec
        [~,l,w] = findpeaks(Zs(:,c),'MinPeakHeight',Opt.DSThrSD,'MinPeakDistance',round(0.1*fsL),'WidthReference','halfheight');
        wMs = w/fsL*1000;
        l = l(wMs >= Opt.DSWidthMs(1) & wMs <= Opt.DSWidthMs(2) & ~nearRip(l));
        dsPk{s,c} = l; nDS(c) = nDS(c) + numel(l);
    end
end
hwD = round(0.1*fsL); tD = (-hwD:hwD)/fsL;
for sh = unique(Shank(:))'
    ch = find(Shank == sh);
    [nRef,k] = max(nDS(ch)); ref = ch(k);
    if nRef < Opt.DSMinEvents, continue; end
    Sum = zeros(numel(ch),numel(tD)); n = 0;
    for s = find(Quiet)'
        for t0 = dsPk{s,ref}(:)'
            if t0 - hwD < 1 || t0 + hwD > size(Ylow{s},1), continue; end
            Sum = Sum + double(Ylow{s}(t0-hwD:t0+hwD,ch))'; n = n + 1;
        end
    end
    C = ComputeCSD(Sum/n,Shank(ch),Depth(ch),tD,struct('BaseWin',[-Inf -0.05]));
    if isempty(C), continue; end
    m = mean(C.CSD(:,abs(tD) <= 0.008),2)/max(C.Noise,eps);
    rd = Depth(ref);
    inSrc = abs(C.Depth - rd) <= 60; inSnk = C.Depth >= rd + 60 & C.Depth <= rd + 350;
    C.RefDepth = rd; C.nEvents = n; C.SourceZ = NaN; C.SinkZ = NaN; C.SourceDepth = NaN; C.SinkDepth = NaN;
    if any(inSrc & ~isnan(m)), [C.SourceZ,i] = max(m.*(inSrc./inSrc)); C.SourceDepth = C.Depth(i); end
    if any(inSnk & ~isnan(m)), [C.SinkZ,i] = min(m.*(inSnk./inSnk)); C.SinkDepth = C.Depth(i); end
    C.Valid = C.SourceZ >= Opt.DSCSDThrSD && C.SinkZ <= -Opt.DSCSDThrSD;
    if isempty(DentateCSD), DentateCSD = C; else, DentateCSD(end+1) = C; end %#ok<AGROW>
end

%% Theta-trough-triggered CSD per shank (channel with the most theta power; cycles with above-median
% theta amplitude). In the hippocampus theta has a laminar sink/source pattern (largest near the
% hippocampal fissure, phase reversal across CA1); for the curator, not fitted
ThetaCSD = [];
fsL = fsOut/dsS; [bT,aT] = butter(2,[6 10]/(fsL/2),'bandpass');
thetaPow = mean(Pxx(f >= 6 & f <= 10,:),1)';
for sh = unique(Shank(:))'
    ch = find(Shank == sh);
    [~,k] = max(thetaPow(ch)); ref = ch(k);
    Sum = zeros(numel(ch),numel(tS)); n = 0;
    for s = find(Active)'
        h = hilbert(filtfilt(bT,aT,double(Ylow{s}(:,ref))));
        ph = angle(h); amp = abs(h);
        tr = find(diff(ph) < -pi) + 1; % phase wraps from +pi to -pi: trough
        tr = tr(amp(tr) > median(amp) & tr - hw >= 1 & tr + hw <= size(Ylow{s},1));
        for t0 = tr(:)'
            Sum = Sum + double(Ylow{s}(t0-hw:t0+hw,ch))'; n = n + 1;
        end
    end
    if n < 50, continue; end
    C = ComputeCSD(Sum/n,Shank(ch),Depth(ch),tS,struct('BaseWin',[-Inf Inf])); % noise: overall SD (theta is ongoing)
    if isempty(C), continue; end
    m = mean(C.CSD(:,abs(tS) <= 0.01),2)/max(C.Noise,eps);
    [C.SinkZ,i] = min(m); C.SinkDepth = C.Depth(i);
    [C.SourceZ,i] = max(m); C.SourceDepth = C.Depth(i);
    C.RefDepth = Depth(ref); C.nEvents = n;
    if isempty(ThetaCSD), ThetaCSD = C; else, ThetaCSD(end+1) = C; end %#ok<AGROW>
end

%% Slow-wave CSD per shank (quiet segments). No assumption on where cortex is on the probe: the CSD
% removes volume conduction, so the depth with the most delta-band (1-4 Hz) CSD power is a local
% generator. Excluding the hippocampus only where the data show it (theta CSD >= Opt.SlowExclThetaSD,
% or >= 10% of its maximum if that is hippocampal-strong (>= Opt.SlowHippoThetaSD); theta CSD is only
% generated there; and +-300 um around a confirmed CA1 ripple source), the depth with the most delta
% CSD power is the reference; the
% CSD is averaged around its large sinks (troughs of its delta CSD, amplitude above the 75th
% percentile). In the neocortex the slow wave (up state) comes with a sink in the deep layers
% (layer 5) and a source above: Valid if that sink/source pair is clear (>= Opt.SlowCSDThrSD).
SlowCSD = [];
[bD,aD] = butter(2,[1 4]/(fsL/2),'bandpass');
hwS = round(0.3*fsL); tSW = (-hwS:hwS)/fsL;
for sh = unique(Shank(:))'
    ch = find(Shank == sh);
    % delta-band CSD time series of each quiet segment
    Cs = cell(numel(starts),1); pw = 0; grid_ = [];
    for s = find(Quiet)'
        X = filtfilt(bD,aD,double(Ylow{s}(:,ch)))';
        Cq = ComputeCSD(X,Shank(ch),Depth(ch),(1:size(X,2))/fsL,struct('BaseWin',[-Inf Inf]));
        if isempty(Cq), continue; end
        Cs{s} = Cq.CSD; grid_ = Cq.Depth; pw = pw + mean(Cq.CSD.^2,2,'omitnan');
    end
    if isempty(grid_), continue; end
    tc = []; if ~isempty(ThetaCSD), tc = ThetaCSD([ThetaCSD.Shank] == sh); end
    if ~isempty(tc) % exclude where the theta CSD shows hippocampus (relative to its maximum, if hippocampal-strong)
        rmsT = interp1(tc(1).Depth,sqrt(mean((tc(1).CSD/max(tc(1).Noise,eps)).^2,2,'omitnan')),grid_,'nearest',0);
        if max(rmsT) >= Opt.SlowHippoThetaSD
            pw(rmsT >= max(Opt.SlowExclThetaSD,0.1*max(rmsT))) = NaN;
        else
            pw(rmsT >= Opt.SlowExclThetaSD) = NaN;
        end
    end
    rc = []; if ~isempty(RippleCSD), rc = RippleCSD([RippleCSD.Shank] == sh & [RippleCSD.Valid]); end
    if ~isempty(rc), pw(abs(grid_ - rc(1).SourceDepth) <= 300) = NaN; end % around confirmed CA1 (ripples arise only there)
    [~,iRef] = max(pw);
    Sum = zeros(numel(grid_),numel(tSW)); n = 0;
    for s = find(Quiet)'
        if isempty(Cs{s}), continue; end
        h = hilbert(Cs{s}(iRef,:)); ph = angle(h); amp = abs(h);
        tr = find(diff(ph) < -pi) + 1; % troughs of the delta CSD: sinks at the reference depth
        tr = tr(amp(tr) > prctile(amp,75) & tr - hwS >= 1 & tr + hwS <= size(Cs{s},2));
        for t0 = tr(:)'
            Sum = Sum + Cs{s}(:,t0-hwS:t0+hwS); n = n + 1;
        end
    end
    if n < 30, continue; end
    A = Sum/n; nz = median(std(A,0,2,'omitnan'),'omitnan'); % CSD (already per grid depth)
    m = mean(A(:,abs(tSW) <= 0.025),2)/max(nz,eps);
    m(isnan(pw)) = NaN; % the same exclusions (hippocampus, probe ends) for the sink
    [sinkZ,i] = min(m); sinkD = grid_(i);
    near = abs(grid_ - sinkD) <= 400; % its dipole partner
    [srcZ,j] = max(m.*(near./near));
    C = struct('Shank',sh,'Depth',grid_,'t',tSW,'LFP',[],'CSD',A,'Noise',nz,'SinkZ',sinkZ,'SinkDepth',sinkD, ...
        'SourceZ',srcZ,'SourceDepth',grid_(j),'RefDepth',grid_(iRef),'nEvents',n, ...
        'Valid',-sinkZ >= Opt.SlowCSDThrSD && srcZ >= Opt.SlowCSDThrSD);
    if isempty(SlowCSD), SlowCSD = C; else, SlowCSD(end+1) = C; end %#ok<AGROW>
end
Pxx = Pxx/nP; NbCorr = NbCorr/numel(starts);
Bands = {'delta 1-4 Hz',[1 4]; 'theta 6-10 Hz',[6 10]; 'gamma 30-80 Hz',[30 80]; 'ripple band 150-250 Hz',[150 250]};
Values = nan(nRec,size(Bands,1)+2);
for b = 1:size(Bands,1)
    in = f >= Bands{b,2}(1) & f <= Bands{b,2}(2);
    Values(:,b) = log10(mean(Pxx(in,:),1))';
end
Values(:,end-1) = nEvents/max(quietSec,eps)*60; % during quiet segments
Values(:,end) = NbCorr;
Values(:,end+1) = nDS/max(quietSec,eps)*60;
LFP = struct('Shank',Shank,'Depth',Depth,'Values',Values,'Names',{[Bands(:,1)' {'ripple events per min','neighbour correlation','dentate spikes per min'}]}, ...
    'File',rawFile,'Seconds',totSec,'RippleCSD',RippleCSD,'ThetaCSD',ThetaCSD,'DentateCSD',DentateCSD,'SlowCSD',SlowCSD,'State',State);
end
