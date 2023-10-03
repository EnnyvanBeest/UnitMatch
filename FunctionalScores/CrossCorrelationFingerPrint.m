function [FingerprintRAll,SigMask,AllSessionCorrelationsFingerprints] = CrossCorrelationFingerPrint(sessionCorrelationsAll,Pairs,Unit2Take,recsesGood,plt)
    %% This function will compute the cross-correlation fingerprint.

    if nargin<5
        plt = 1;
    end

    %% Parameters
    nclus = numel(Unit2Take);
    ndays = numel(sessionCorrelationsAll);
    RecOpt = unique(recsesGood);
    SessionSwitch = [1 1+cumsum(cell2mat(cellfun(@(x) size(x.fold1,1), sessionCorrelationsAll, 'uni', 0)))];
    AllSessionCorrelationsFingerprints = cell(ndays,ndays); % Used for plotting only
    FingerprintR_plt = cell(ndays,ndays);
    FingerprintRAll = nan(nclus,nclus);
    SigMask = zeros(nclus,nclus);

    %% Get average cross-correlation matrix
    SessionCorrelations_perDay = cell(1,ndays);
    for did1 = 1:ndays
        % Get average for non x-validated variables (z-transform mean)
        SessionCorrelations_perDay{did1} = tanh(nanmean(atanh(cat(3,sessionCorrelationsAll{did1}.fold1,sessionCorrelationsAll{did1}.fold2)),3));
    end

    %% Computes the fingerprint correlations and rank across days
    for did1 = 1:ndays
        for did2 = 1:ndays
            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1)-1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1)-1;
            if did2==did1
                % Within day, compute correlations across the 2 folds

                % Get the correlation of the fingerprints
                FingerprintR = corr(sessionCorrelationsAll{did1}.fold1,sessionCorrelationsAll{did1}.fold2,'rows','pairwise'); % cross-validated

                % Save correlation matrix
                AllSessionCorrelationsFingerprints{did1,did2} = cat(2,sessionCorrelationsAll{did1}.fold1,sessionCorrelationsAll{did1}.fold2)'; % for plotting save both folds

            else
                % Across days, compute correlations across days

                dayopt = [did1,did2];
                SessionCorrelation_Pair = cell(1,2);
                
                % Find the pairs
                for did = 1:length(dayopt)
                    % We need a group of units that is likely to be a pair across at least two days
                    if did==1
                        pairidx = recsesGood(Pairs(:,1)) == RecOpt(dayopt(did)) & recsesGood(Pairs(:,2))==RecOpt(dayopt(did+1));
                        PairsTmp = Pairs(pairidx,:);
                        % Only use every 'unit' once --> take the highest scoring matches
                        [~,id1,~]=unique(PairsTmp(:,1),'stable');
                        PairsTmp = PairsTmp(id1,:);
                        [~,id1,~]=unique(PairsTmp(:,2),'stable');
                        PairsTmp = PairsTmp(id1,:);
                        Unit2TakeIdx = PairsTmp(:,1); % Only take each unit once
                    else
                        Unit2TakeIdx = [];
                    end
                    if did==2
                        pairidx = recsesGood(Pairs(:,2)) == RecOpt(dayopt(did)) & recsesGood(Pairs(:,1))==RecOpt(dayopt(did-1));
                        PairsTmp = Pairs(pairidx,:);
                        % Only use every 'unit' once --> take the highest scoring matches
                        [~,id1,~]=unique(PairsTmp(:,1),'stable');
                        PairsTmp = PairsTmp(id1,:);
                        [~,id1,~]=unique(PairsTmp(:,2),'stable');
                        PairsTmp = PairsTmp(id1,:);
                        Unit2TakeIdx = [Unit2TakeIdx; PairsTmp(:,2)];
                    end
                  
                    Unit2TakeIdxAll = find(recsesGood == RecOpt(dayopt(did)));
    
                    % Extract the part of the correlation matrix with these
                    % pairs
                    sortIdx = cell2mat(arrayfun(@(x) find(Unit2Take(Unit2TakeIdxAll) == x), Unit2Take(Unit2TakeIdx), 'uni', 0));
                    SessionCorrelation_Pair{did} = SessionCorrelations_perDay{dayopt(did)}(sortIdx,:);

                    % Remove =1 for the same unit (per def. 1)
                    for id = 1:length(Unit2TakeIdx)
                        SessionCorrelation_Pair{did}(id,ismember(Unit2TakeIdxAll,Unit2TakeIdx(id))) = nan;
                    end
    
                    % Normalize correlations to compare across recordings
%                     SessionCorrelation_Pair{did} = (SessionCorrelation_Pair{did}-nanmedian(SessionCorrelation_Pair{did}(:)))./ ...
%                         (quantile(SessionCorrelation_Pair{did}(:),0.95)-quantile(SessionCorrelation_Pair{did}(:),0.05));
                end

                % Get the correlation of the fingerprints
                if isempty(SessionCorrelation_Pair{1}) || isempty(SessionCorrelation_Pair{2})
                    FingerprintR = zeros(length(clusIdxD1All),length(clusIdxD2All));
                else
                    FingerprintR = corr(SessionCorrelation_Pair{1},SessionCorrelation_Pair{2},'rows','pairwise');
                end

                AllSessionCorrelationsFingerprints{did1,did2} = cat(2,SessionCorrelation_Pair{:})';
            end


            % Save Fingerprint correlations
            FingerprintRAll(clusIdxD1All,clusIdxD2All) = FingerprintR;
        end
    end

    %% Plots
    if plt

        %% Plot the correlation matrices

        figure('name','Cross-correlation Fingerprints');
        nrows = (ndays*ndays);
        rowcount = 1;

        tic
        for did1 = 1:ndays
            for did2 = 1:ndays
                if did2==did1
                    subplot(nrows,2,(rowcount-1)*2+1)
                    imagesc(sessionCorrelationsAll{did1}.fold1')
                    colormap(flipud(gray))
                    xlabel('Candidate Units to be matched')
                    ylabel(['Within day ' num2str(did1)])
                    title('First half of recording')
                    makepretty

                    subplot(nrows,2,rowcount*2)
                    imagesc(sessionCorrelationsAll{did1}.fold2')
                    colormap(flipud(gray))
                    xlabel('Candidate Units to be matched')
                    ylabel(['Within day ' num2str(did1)])
                    title('Second half of recording')
                    makepretty
                else
                    clusIdxD{1} = 1:diff(SessionSwitch(did1+(0:1)));
                    clusIdxD{2} = diff(SessionSwitch(did1+(0:1)))+(1:diff(SessionSwitch(did2+(0:1))));
                    for did = 1:2
                        subplot(nrows,2,(rowcount-1)*2+did)
                        imagesc(AllSessionCorrelationsFingerprints{did1,did2}(clusIdxD{did},:)')
                        colormap(flipud(gray))
                        xlabel('Candidate Units to be matched')
                        ylabel(['Across days ' num2str(dayopt(did))])
                        title(['Recording ' num2str(dayopt(did))])
                        makepretty
                    end
                end
                rowcount = rowcount+1;
            end
        end

        %% Plot the correlation of the fingerprints

        figure('name','Fingerprint correlations')
        for did1 = 1:ndays
            for did2 = 1:ndays
                subplot(ndays,ndays,(did1-1)*ndays+did2)
                clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1)-1;
                clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1)-1;

                imagesc(FingerprintRAll(clusIdxD1All,clusIdxD2All))
                hold on
                colormap(flipud(gray))
                xlabel('All units across both days')
                ylabel('All units across both days')
                title(['R' num2str(did2)])
                ylabel(['R' num2str(did1)])
                set(gca,'XTickLabel','','YTickLabel','')
                makepretty
            end
        end
    end
end

%% Back-up
% 
% %% Plot the rank
% 
% figure('name','RankScore')
% subplot(1,2,1)
% imagesc(RankScoreAll==1)
% hold on
% arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
% arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
% colormap(flipud(gray))
% xlabel('All units across both days')
% ylabel('All units across both days')
% title('RankScore = 1')
% makepretty
% subplot(1,2,2)
% imagesc(SigMask)
% hold on
% arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
% arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
% colormap(flipud(gray))
% xlabel('All units across both days')
% ylabel('All units across both days')
% title('correlations>99th percentile of distribution')
% makepretty