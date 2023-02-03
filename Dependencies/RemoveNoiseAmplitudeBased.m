function sp = RemoveNoiseAmplitudeBased(sp)

disp('Removing noise spikes...')
figure
subplot(1,2,1)
scatter(sp.st,sp.spikeDepths,4,[0 0 0],'filled')
title('Spikes with noise')
hold on

UniqClus = unique(sp.spikeTemplates);
rmidx = [];
PercRemoved = nan(1,length(UniqClus));
for uid = 1:length(UniqClus)
    tmpfind = find(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid))>quantile(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid)),0.75)+1.5*(quantile(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid)),0.75)-quantile(sp.spikeAmps(sp.spikeTemplates==UniqClus(uid)),0.25)));
    tmpidx = find(sp.spikeTemplates==UniqClus(uid));
    PercRemoved(uid) = length(tmpfind)./length(tmpidx)*100;
    rmidx = [rmidx; tmpidx(tmpfind)];
end
scatter(sp.st(rmidx),sp.spikeDepths(rmidx),4,[1 0 0],'filled')

drawnow
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
scatter(sp.st,sp.spikeDepths,4,[0 0 0],'filled')
title('Spikes after noise removal')
hold on


return