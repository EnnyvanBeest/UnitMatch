function sp = RemoveNoiseAmplitudeBased(sp)
RemoveIndividualNoise = 0; % Individual unit noise --> outlier amplitudes
RemoveCommonNoise = 1; % Noise shared across channels (e.g. lick respones)

figure
subplot(1,2,1)
scatter(sp.st,sp.spikeDepths,1,[0 0 0],'filled')
title('Spikes with noise')
hold on

if RemoveCommonNoise
    disp('Removing common noise...')
    % First remove general noise
    binsz = 0.001;
    edges = [min(sp.st)-binsz/2:binsz:max(sp.st)+binsz/2];
    timevec = [min(sp.st):binsz:max(sp.st)];

    depthstep = 25; %um
    depthedges = [min(sp.spikeDepths)-depthstep/2:depthstep:max(sp.spikeDepths)+depthstep/2];
    tmpact = arrayfun(@(X) histcounts(sp.st(sp.spikeDepths>depthedges(X)&sp.spikeDepths<depthedges(X+1)),edges),1:length(depthedges)-1,'UniformOutput',0);
    tmpact = cat(1,tmpact{:});
    tmpact = (tmpact-nanmean(tmpact,2))./nanstd(tmpact,[],2); %Z-score
    tpidx = find(sum(tmpact>3,1)>0.5*length(depthedges)-1);

    % Remove these timepoints
    rmidx = arrayfun(@(X) find(sp.st>timevec(X)-binsz/2&sp.st<timevec(X)+binsz/2),tpidx,'UniformOutput',0);
    rmidx = cat(1,rmidx{:});
end

if RemoveIndividualNoise
    % Remove Individual Noise
    disp('Removing individual unit noise')
    % Now remove noisy spikes of individual units
    UniqClus = unique(sp.spikeTemplates);
    PercRemoved = nan(1,length(UniqClus));
    for uid = 1:length(UniqClus)
        tmpfind = find(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid))>quantile(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid)),0.75)+1.5*(quantile(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid)),0.75)-quantile(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid)),0.25)));
        tmpidx = find(sp.spikeTemplates==UniqClus(uid));
        PercRemoved(uid) = length(tmpfind)./length(tmpidx)*100;
        rmidx = [rmidx; tmpidx(tmpfind)];
    end

end
scatter(sp.st(rmidx),sp.spikeDepths(rmidx),1,[1 0 0],'filled')

nori = length(sp.st);
fields = fieldnames(sp);
for fid = 1:length(fields)
    eval(['tmp = sp. ' fields{fid} ';'])
    if any(size(tmp) == nori)
        tmp(rmidx,:,:)=[];
        eval(['sp.' fields{fid} '=tmp;'])
    end
end

subplot(1,2,2)
scatter(sp.st,sp.spikeDepths,1,[0 0 0],'filled')
title('Spikes after noise removal')
hold on
drawnow

sp.noisespkes = rmidx;
return