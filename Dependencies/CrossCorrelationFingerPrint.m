timevec = floor(min(sp.st)):binsz:ceil(max(sp.st));
edges = floor(min(sp.st))-binsz/2:binsz:ceil(max(sp.st))+binsz/2;
disp('Calculate activity correlations')

%% Correlations per session
Unit2Take = AllClusterIDs(Good_Idx);
fig1 = figure('name','Cross-correlation Fingerprints');
clear AllSessionCorrelations
if ndays>1
    for did1 = 1:ndays-1
        for did2 = 2:ndays
            if did2<=did1
                continue
            end
            dayopt = [did1,did2];
            for did = 1:length(dayopt)
                PairsTmp = Pairs;
                % We need a group of units that is likely to be a pair across at least two days
                if did==1
                    pairidx = find(recsesAll(Good_Idx(Pairs(:,1))) == dayopt(did) & recsesAll(Good_Idx(Pairs(:,2)))==dayopt(did+1));
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
                if did==2
                    pairidx = find(recsesAll(Good_Idx(Pairs(:,2))) == dayopt(did) & recsesAll(Good_Idx(Pairs(:,1)))==dayopt(did-1));
                    PairsTmp = Pairs(pairidx,:);
                    % Only use every 'unit' once --> take the highest scoring matches
                    [val,id1,id2]=unique(PairsTmp(:,1),'stable');
                    PairsTmp = PairsTmp(id1,:);
                    [val,id1,id2]=unique(PairsTmp(:,2),'stable');
                    PairsTmp = PairsTmp(id1,:);
                    Unit2TakeIdx = [Unit2TakeIdx; PairsTmp(:,2)];
                end
                % Correlation on this day
                srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == dayopt(did)),edges),Unit2Take(Unit2TakeIdx),'UniformOutput',0);
                srMatches = cat(1,srMatches{:});

                % All Units on this day
                Unit2TakeIdxAll = find(recsesAll(Good_Idx) == dayopt(did));
                srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == dayopt(did)),edges),Unit2Take(Unit2TakeIdxAll),'UniformOutput',0);
                srAll = cat(1,srAll{:});
                SessionCorrelation_Pair = corr(srMatches(:,1:floor(size(srMatches,2)./2))',srAll(:,1:floor(size(srMatches,2)./2))');

                %     % Remove =1 for the same unit (per def. 1)
                for id = 1:length(Unit2TakeIdx)
                    SessionCorrelation_Pair(id,find(ismember(Unit2TakeIdxAll,Unit2TakeIdx(id)))) = nan;
                end

                % Normalize correlations to compare across recordings
                %         SessionCorrelation_Pair = (SessionCorrelation_Pair-nanmedian(SessionCorrelation_Pair(:)))./(quantile(SessionCorrelation_Pair(:),0.95)-quantile(SessionCorrelation_Pair(:),0.05));
                subplot(ndays-1,2,did)
                imagesc(SessionCorrelation_Pair')
                colormap(flipud(gray))
                xlabel('Candidate Units to be matched')
                ylabel('All units')
                title(['Recording ' num2str(dayopt(did))])
                makepretty

                % Add all together
                if did == 1
                    SessionCorrelations = SessionCorrelation_Pair';
                else
                    SessionCorrelations = cat(1,SessionCorrelations,SessionCorrelation_Pair');
                end
            end
            AllSessionCorrelations{did1,did2} = SessionCorrelations;
        end
    end
else
    % All Units on this day
    Unit2TakeIdxAll = find(recsesAll(Good_Idx) == 1);
    srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == 1),edges),Unit2Take(Unit2TakeIdxAll),'UniformOutput',0);
    srAll = cat(1,srAll{:});

    % Find cross-correlation in first and second half of session
    SessionCorrelations = corr(srAll(:,1:floor(size(srAll,2)./2))',srAll(:,1:floor(size(srAll,2)./2))')';
    SessionCorrelations2 = corr(srAll(:,floor(size(srAll,2)./2)+1:floor(size(srAll,2)./2)*2)',srAll(:,floor(size(srAll,2)./2)+1:floor(size(srAll,2)./2)*2)')';

    %     % Remove =1 for the same unit (per def. 1)
    for id = 1:length(Unit2TakeIdxAll)
        SessionCorrelations(id,find(ismember(Unit2TakeIdxAll,Unit2TakeIdxAll(id)))) = nan;
        SessionCorrelations2(id,find(ismember(Unit2TakeIdxAll,Unit2TakeIdxAll(id)))) = nan;
    end

    subplot(1,2,1)
    imagesc(SessionCorrelations')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title(['First half of recording'])
    makepretty

    subplot(1,2,2)
    imagesc(SessionCorrelations2')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title(['Second half of recording'])
    makepretty


    SessionCorrelations = nanmean(cat(3,SessionCorrelations',SessionCorrelations2'),3);
    AllSessionCorrelations{1} = SessionCorrelations;
end

%%
ncellsperrecording = diff(SessionSwitch);

clear FingerPrintAll
figure('name','Fingerprint correlations')
for did1 = 1:ndays
    for did2 = 1:ndays
        if did2<=did1 && ndays~=1
            continue
        end
        SessionCorrelations = AllSessionCorrelations{did1,did2};
        rmidx = find(sum(isnan(SessionCorrelations),2)==size(SessionCorrelations,2));
        nclustmp = size(SessionCorrelations,1);
%         SessionCorrelations(rmidx,:)=[];
        try
            notrmdixvec = 1:ncellsperrecording(did1)+ncellsperrecording(did2);
            notrmdixvec(rmidx)=[];
        catch
            notrmdixvec = SessionSwitch(did1):SessionSwitch(did1+1)-1;
            notrmdixvec(rmidx)=[];

        end
        % Correlate 'fingerprints'
        FingerprintR = arrayfun(@(X) cell2mat(arrayfun(@(Y) corr(SessionCorrelations(X,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))',SessionCorrelations(Y,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))'),notrmdixvec,'UniformOutput',0)),notrmdixvec,'UniformOutput',0);
        FingerprintR = cat(1,FingerprintR{:});

        % If one value was only nans; put it back in and replace original
        % FingerprintR
        if any(rmidx)
            Fingerprinttmp = nan(nclustmp,nclustmp);
            Fingerprinttmp(1:rmidx-1,1:rmidx-1)=FingerprintR(1:rmidx-1,1:rmidx-1);
            Fingerprinttmp(rmidx,:)=nan;
            Fingerprinttmp(:,rmidx)=nan;
            Fingerprinttmp(rmidx+1:end,notrmdixvec)=FingerprintR(rmidx:end,:);
            Fingerprinttmp(notrmdixvec,rmidx+1:end)=FingerprintR(:,rmidx:end);
            FingerprintR = Fingerprinttmp;
            clear Fingerprinttmp
        end

        if ndays>1
            subplot(ndays-1,ndays-1,(did1-1)*(ndays-1)+did2-1)
        end
        imagesc(FingerprintR)
        hold on
        arrayfun(@(X) line([SessionSwitch(X)-SessionSwitch(X-1) SessionSwitch(X)-SessionSwitch(X-1)],get(gca,'ylim'),'color',[1 0 0]),did1+1,'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X)-SessionSwitch(X-1) SessionSwitch(X)-SessionSwitch(X-1)],'color',[1 0 0]),did1+1,'Uni',0)
        colormap(flipud(gray))
%         xlabel('All units across both days')
%         ylabel('All units across both days')
        title(['R' num2str(did2)])
        ylabel(['R' num2str(did1)])
        set(gca,'XTickLabel','','YTickLabel','')

        makepretty
        drawnow

        FingerPrintAll{did1,did2} = FingerprintR;
    end
end
%%
FingerprintRAll = nan(nclus,nclus);
SigMask = zeros(nclus,nclus);
RankScoreAll = nan(size(SigMask));
for pid=1:nclus
    for pid2 = 1:nclus
        if pid2<pid
            continue
        end
        did1 = (recsesGood(pid));
        did2 = (recsesGood(pid2));
        addthis3=-SessionSwitch(did1)+1;
        addthis4 = 0;
        if did2==did1 && did2~=ndays
            did2=did2+1;
        elseif did2==did1 && did2==ndays
            try
                addthis3= -SessionSwitch(did1)+1+ncellsperrecording(did1-1); %for pid we need to subtract when there's mltiple sessions
                addthis4 = ncellsperrecording(did1-1); % For the columns (pid2)
                did1 = did1-1;           % Subtract a day, because we didn't redo this up there
            catch
                addthis3= -SessionSwitch(did1)+1; %for pid we need to subtract when there's mltiple sessions
                addthis4 = 0; % For the columns (pid2)
            end
        end
      
        FingerprintR = FingerPrintAll{did1,did2};

        % Reset days
        did2 = (recsesGood(pid2));
        did1 = (recsesGood(pid));

        if did1==did2
            tmp1 = FingerprintR(pid+addthis3,[1:ncellsperrecording(did1)]+addthis4);
            addthis=SessionSwitch(did1)-1;
            addthis2 = addthis3;
        else
            tmp1 = FingerprintR(pid+addthis3,ncellsperrecording(did1)+1:end);
            addthis=SessionSwitch(did2)-1;
            addthis2 = -addthis+ncellsperrecording(did1);
        end
        tmp1(isnan(tmp1))=0;
        [val,ranktmp] = sort(tmp1,'descend');

        tmp1(pid2-addthis)=[];

        if FingerprintR(pid+addthis3,pid2+addthis2)>quantile(tmp1,0.99)
            SigMask(pid,pid2)=1;
        end
        FingerprintRAll(pid,pid2) = FingerprintR(pid+addthis3,pid2+addthis2);
        RankScoreAll(pid,pid2) = find(ranktmp==pid2-addthis);
    end
end

% MIRROR
for uid=1:nclus
    for uid2 = 1:nclus
        if uid2<uid
            if ~isnan(RankScoreAll(uid,uid2))
                keyboard
            end
            SigMask(uid,uid2)=SigMask(uid2,uid);
            RankScoreAll(uid,uid2)=RankScoreAll(uid2,uid);
            FingerprintRAll(uid,uid2) = FingerprintRAll(uid2,uid);
        end
    end
end

FingerprintR = FingerprintRAll;
clear FingerprintRAll

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
title('correlations>99th percentile of distribution')
makepretty