timevec = floor(min(sp.st)):binsz:ceil(max(sp.st));
edges = floor(min(sp.st))-binsz/2:binsz:ceil(max(sp.st))+binsz/2;
disp('Calculate activity correlations')

%% Correlations per session
if ndays>1
    Unit2Take = AllClusterIDs(Good_Idx);
    fig1 = figure('name','Cross-correlation Fingerprints');
    for did = 1:ndays
        PairsTmp = Pairs;
        % We need a group of units that is likely to be a pair across at least two days
        if did<ndays
            pairidx = find(recsesAll(Good_Idx(Pairs(:,1))) == did & recsesAll(Good_Idx(Pairs(:,2)))==did+1);
            PairsTmp = Pairs(pairidx,:);
            % Only use every 'unit' once --> take the highest scoring matches
            [val,id1,id2]=unique(PairsTmp(:,1),'stable');
            PairsTmp = PairsTmp(id1,:);
            [val,id1,id2]=unique(PairsTmp(:,2),'stable');
            PairsTmp = PairsTmp(id1,:);
            Unit2TakeIdx = PairsTmp(:,1); % Only take each unit once
        else
            Unit2TakeIdx = [];
        end
        if did>1
            pairidx = find(recsesAll(Good_Idx(Pairs(:,2))) == did & recsesAll(Good_Idx(Pairs(:,1)))==did-1);
            PairsTmp = Pairs(pairidx,:);
            % Only use every 'unit' once --> take the highest scoring matches
            [val,id1,id2]=unique(PairsTmp(:,1),'stable');
            PairsTmp = PairsTmp(id1,:);
            [val,id1,id2]=unique(PairsTmp(:,2),'stable');
            PairsTmp = PairsTmp(id1,:);
            Unit2TakeIdx = [Unit2TakeIdx; PairsTmp(:,2)];
        end
        % Correlation on this day
        srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == did),edges),Unit2Take(Unit2TakeIdx),'UniformOutput',0);
        srMatches = cat(1,srMatches{:});

        % All Units on this day
        Unit2TakeIdxAll = find(recsesAll(Good_Idx) == did);
        srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == did),edges),Unit2Take(Unit2TakeIdxAll),'UniformOutput',0);
        srAll = cat(1,srAll{:});
        SessionCorrelation_Pair = corr(srMatches(:,1:floor(size(srMatches,2)./2))',srAll(:,1:floor(size(srMatches,2)./2))');

        %     % Remove =1 for the same unit (per def. 1) - no longer true, use
        %     first and second half of recording
        for id = 1:length(Unit2TakeIdx)
            SessionCorrelation_Pair(id,find(ismember(Unit2TakeIdxAll,Unit2TakeIdx(id)))) = nan;
        end

        % Normalize correlations to compare across recordings
%         SessionCorrelation_Pair = (SessionCorrelation_Pair-nanmedian(SessionCorrelation_Pair(:)))./(quantile(SessionCorrelation_Pair(:),0.95)-quantile(SessionCorrelation_Pair(:),0.05));
        subplot(1,ndays,did)
        imagesc(SessionCorrelation_Pair')
        colormap(flipud(gray))
        xlabel('Candidate Units to be matched')
        ylabel('All units')
        title(['Recording ' num2str(did)])
        makepretty

        % Add all together
        if did == 1
            SessionCorrelations = SessionCorrelation_Pair';
        else
            SessionCorrelations = cat(1,SessionCorrelations,SessionCorrelation_Pair');
        end
    end
else
    PairsTmp = Pairs;
    % All Units on this day
    Unit2TakeIdxAll = find(recsesAll(Good_Idx) == 1);
    srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == 1),edges),Unit2Take(Unit2TakeIdxAll),'UniformOutput',0);
    srAll = cat(1,srAll{:});

    SessionCorrelations = corr(srAll(:,1:floor(size(srAll,2)./2))',srAll(:,floor(size(srAll,2)./2)+1:floor(size(srAll,2)./2)*2)')';
    keyboard %Check this
    subplot(1,ndays,did)
    imagesc(SessionCorrelations')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title(['Recording ' num2str(did)])
    makepretty
end

%%
rmidx = find(sum(isnan(SessionCorrelations),2)==size(SessionCorrelations,2));
SessionCorrelations(rmidx,:)=[];
nclustmp = nclus-length(rmidx);
notrmdixvec = 1:nclus;
notrmdixvec(rmidx)=[];
% Correlate 'fingerprints'
FingerprintR = arrayfun(@(X) cell2mat(arrayfun(@(Y) corr(SessionCorrelations(X,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))',SessionCorrelations(Y,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))'),1:nclustmp,'UniformOutput',0)),1:nclustmp,'UniformOutput',0);
FingerprintR = cat(1,FingerprintR{:});

% If one value was only nans; put it back in and replace original
% FingerprintR
if any(rmidx)
    Fingerprinttmp = nan(nclus,nclus);
    Fingerprinttmp(1:rmidx-1,1:rmidx-1)=FingerprintR(1:rmidx-1,1:rmidx-1);
    Fingerprinttmp(rmidx,:)=nan;
    Fingerprinttmp(:,rmidx)=nan;
    Fingerprinttmp(rmidx+1:end,notrmdixvec)=FingerprintR(rmidx:end,:);
    Fingerprinttmp(notrmdixvec,rmidx+1:end)=FingerprintR(:,rmidx:end);
    FingerprintR = Fingerprinttmp;
    clear Fingerprinttmp
end

figure('name','Fingerprint correlations')
imagesc(FingerprintR)
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
xlabel('All units across both days')
ylabel('All units across both days')
title('Correlation Fingerprint')
makepretty

%%
SigMask = zeros(nclus,nclus);
RankScoreAll = nan(size(SigMask));
for pid=1:nclus
    for pid2 = 1:nclus
        tmp1 = FingerprintR(pid,SessionSwitch(recsesGood(pid2)):SessionSwitch(recsesGood(pid2)+1)-1);
        addthis=SessionSwitch(recsesGood(pid2))-1;

        [val,ranktmp] = sort(tmp1,'descend');

        tmp1(pid2-addthis)=[];

        if FingerprintR(pid,pid2)>nanmean(tmp1)+2*nanstd(tmp1)
            SigMask(pid,pid2)=1;
        end
        RankScoreAll(pid,pid2) = find(ranktmp==pid2-addthis);
    end
end

figure('name','RankScore')
subplot(1,2,1)
imagesc(RankScoreAll==1)
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
xlabel('All units across both days')
ylabel('All units across both days')
title('RankScore = 1')
makepretty
subplot(1,2,2)
imagesc(SigMask)
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
xlabel('All units across both days')
ylabel('All units across both days')
title('correlations>95th percentile of distribution')
makepretty