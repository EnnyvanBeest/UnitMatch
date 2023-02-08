timevec = floor(min(sp.st)):binsz:ceil(max(sp.st));
edges = floor(min(sp.st))-binsz/2:binsz:ceil(max(sp.st))+binsz/2;

%% Correlations per session
Unit2Take = AllClusterIDs(Good_Idx);
fig1 = figure('name','Cross-correlation Fingerprints');
nrows = (ndays*ndays);
rowcount = 1;
AllSessionCorrelations = cell(ndays,ndays);

% srAllDays = cell(1,ndays); 
% for did = 1:ndays
%     Unit2TakeIdxAll = find(recsesAll(Good_Idx) == did);
%     srAllDays{did} = nan(numel(Unit2TakeIdxAll),numel(edges)-1);
%     for uid = 1:numel(Unit2TakeIdxAll)
%         srAllDays{did}(uid,:) =  histcounts(sp.st(sp.spikeTemplates == Unit2Take(Unit2TakeIdxAll(uid)) & sp.RecSes == did),edges);
%     end
% end

for did1 = 1:ndays
    for did2 = did1:ndays
        if did2==did1
            % All Units on this day
            Unit2TakeIdxAll = find(recsesAll(Good_Idx) == did1);
%             srAll = srAllDays{did1};
            srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == did1),edges),Unit2Take(Unit2TakeIdxAll),'UniformOutput',0);
            srAll = cat(1,srAll{:});
            srAll(:,sum(srAll==0,1)==size(srAll,1))=[]; % Remove bins with 0 spikes

            % Find cross-correlation in first and second half of session
            SessionCorrelations = corr(srAll(:,1:floor(size(srAll,2)./2))',srAll(:,1:floor(size(srAll,2)./2))')';
            SessionCorrelations2 = corr(srAll(:,floor(size(srAll,2)./2)+1:floor(size(srAll,2)./2)*2)',srAll(:,floor(size(srAll,2)./2)+1:floor(size(srAll,2)./2)*2)')';

            %     % Remove =1 for the same unit (per def. 1)
            SessionCorrelations(logical(eye(size(SessionCorrelations)))) = nan;
            SessionCorrelations2(logical(eye(size(SessionCorrelations)))) = nan;

            subplot(nrows,2,(rowcount-1)*2+1)
            imagesc(SessionCorrelations')
            colormap(flipud(gray))
            xlabel('Candidate Units to be matched')
            ylabel(['Within day ' num2str(did1)])
            title(['First half of recording'])
            makepretty

            subplot(nrows,2,rowcount*2)
            imagesc(SessionCorrelations2')
            colormap(flipud(gray))
            xlabel('Candidate Units to be matched')
            ylabel(['Within day ' num2str(did1)])
            title(['Second half of recording'])
            makepretty


            SessionCorrelations = nanmean(cat(3,SessionCorrelations',SessionCorrelations2'),3);
            AllSessionCorrelations{did1,did2} = SessionCorrelations;

        else
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
                Unit2TakeIdxAll = find(recsesAll(Good_Idx) == dayopt(did));

                % Correlation on this day
                sortIdx = cell2mat(arrayfun(@(x) find(Unit2Take(Unit2TakeIdxAll) == x), Unit2Take(Unit2TakeIdx), 'uni', 0));
%                 srMatches = srAllDays{dayopt(did)}(sortIdx,:); % hacky
                srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == dayopt(did)),edges),Unit2Take(Unit2TakeIdx),'UniformOutput',0);
                srMatches = cat(1,srMatches{:});

                % All Units on this day
%                 srAll = srAllDays{dayopt(did)};
                srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == X & sp.RecSes == dayopt(did)),edges),Unit2Take(Unit2TakeIdxAll),'UniformOutput',0);
                srAll = cat(1,srAll{:});
                SessionCorrelation_Pair = corr(srMatches(:,1:floor(size(srMatches,2)./2))',srAll(:,1:floor(size(srMatches,2)./2))');

                %     % Remove =1 for the same unit (per def. 1)
                for id = 1:length(Unit2TakeIdx)
                    SessionCorrelation_Pair(id,find(ismember(Unit2TakeIdxAll,Unit2TakeIdx(id)))) = nan;
                end

                % Normalize correlations to compare across recordings
                %         SessionCorrelation_Pair = (SessionCorrelation_Pair-nanmedian(SessionCorrelation_Pair(:)))./(quantile(SessionCorrelation_Pair(:),0.95)-quantile(SessionCorrelation_Pair(:),0.05));
                subplot(nrows,2,(rowcount-1)*2+did)
                imagesc(SessionCorrelation_Pair')
                colormap(flipud(gray))
                xlabel('Candidate Units to be matched')
                ylabel(['Across days ' num2str(dayopt(did))])
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
        rowcount=rowcount+1;
    end
end

%%
ncellsperrecording = diff(SessionSwitch);
rowcount = 1;
clear FingerPrintAll
figure('name','Fingerprint correlations')
for did1 = 1:ndays
    for did2 = 1:ndays
        if did2<did1
            continue
        end
        SessionCorrelations = AllSessionCorrelations{did1,did2};
        rmidx = find(sum(isnan(SessionCorrelations),2)==size(SessionCorrelations,2));
        nclustmp = size(SessionCorrelations,1);
        %         SessionCorrelations(rmidx,:)=[];
        try
            notrmdixvec = 1:nclustmp;
            notrmdixvec(rmidx)=[];
        catch
            notrmdixvec = SessionSwitch(did1):SessionSwitch(did1+1)-1;
            notrmdixvec(rmidx)=[];

        end
        % Correlate 'fingerprints'
        try
            x = SessionCorrelations(notrmdixvec,:)';
%             FingerprintR = corr(x,x,'rows','pairwise');
            FingerprintR = arrayfun(@(X) cell2mat(arrayfun(@(Y) corr(SessionCorrelations(X,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))',SessionCorrelations(Y,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))'),notrmdixvec,'UniformOutput',0)),notrmdixvec,'UniformOutput',0);
            FingerprintR = cat(1,FingerprintR{:});
        catch ME
            disp(ME)
            keyboard
        end

        % If one value was only nans; put it back in and replace original
        % FingerprintR
        if any(rmidx)
            Fingerprinttmp = nan(nclustmp,nclustmp);
            Fingerprinttmp(setdiff(1:nclustmp,rmidx),setdiff(1:nclustmp,rmidx)) = FingerprintR;
            FingerprintR = Fingerprinttmp;
            clear Fingerprinttmp
        end

        subplot(ndays,ndays,(did1-1)*ndays+did2)
        imagesc(FingerprintR)
        hold on
        arrayfun(@(X) line([SessionSwitch(X)-SessionSwitch(X-1)+1 SessionSwitch(X)-SessionSwitch(X-1)+1],get(gca,'ylim'),'color',[1 0 0]),did1+1,'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X)-SessionSwitch(X-1)+1 SessionSwitch(X)-SessionSwitch(X-1)+1],'color',[1 0 0]),did1+1,'Uni',0)
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
    
        FingerprintR = FingerPrintAll{did1,did2};
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