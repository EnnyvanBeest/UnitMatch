function [FingerprintRAll,RankScoreAll,SigMask,AllSessionCorrelations] = CrossCorrelationFingerPrint(srAllDays,Pairs,Unit2Take,recsesGood,plt)
    %% This function will compute the cross-correlation fingerprint.

    if nargin<5
        plt = 1;
    end

    %% Parameters

    nclus = numel(Unit2Take);
    ndays = numel(srAllDays);
    SessionSwitch = [1 1+cumsum(cell2mat(cellfun(@(x) size(x,1), srAllDays, 'uni', 0)))];
    AllSessionCorrelations = cell(ndays,ndays); % Used for plotting only
    FingerprintR_plt = cell(ndays,ndays);
    FingerprintRAll = nan(nclus,nclus);
    SigMask = zeros(nclus,nclus);
    RankScoreAll = nan(nclus,nclus);
    
    %% Computes all the cross-correlation matrices
    
    SessionCorrelations_fold1 = cell(1,ndays);
    SessionCorrelations_fold2 = cell(1,ndays);
    SessionCorrelations_perDay = cell(1,ndays);
    for did1 = 1:ndays
        sr = srAllDays{did1};

        % Define folds (two halves)
        idx_fold1 = 1:floor(size(sr,2)./2);
        idx_fold2 = floor(size(sr,2)./2)+1:floor(size(sr,2)./2)*2;

        % Find cross-correlation in first and second half of session
        SessionCorrelations_fold1{did1} = corr(sr(:,idx_fold1)',sr(:,idx_fold1)')';
        SessionCorrelations_fold2{did1} = corr(sr(:,idx_fold2)',sr(:,idx_fold2)')';
        
        % Nan the diagonal
        SessionCorrelations_fold1{did1}(logical(eye(size(SessionCorrelations_fold1{did1})))) = nan;
        SessionCorrelations_fold2{did1}(logical(eye(size(SessionCorrelations_fold2{did1})))) = nan;

        % Get average for non x-validated variables (z-transform mean)
        SessionCorrelations_perDay{did1} = tanh(nanmean(atanh(cat(3,SessionCorrelations_fold1{did1},SessionCorrelations_fold2{did1})),3));
    end

    %% Computes the fingerprint correlations and rank across days

    for did1 = 1:ndays
        for did2 = 1:ndays
            if did2==did1
                % Within day, compute correlations across the 2 folds

                % Get the correlation of the fingerprints
                FingerprintR = corr(SessionCorrelations_fold1{did1},SessionCorrelations_fold2{did1},'rows','pairwise'); % cross-validated

                % Save correlation matrix
                AllSessionCorrelations{did1,did2} = SessionCorrelations_perDay{did1}; % for plotting

            else
                % Across days, compute correlations across days

                dayopt = [did1,did2];
                SessionCorrelation_Pair = cell(1,2);
                
                % Find the pairs
                for did = 1:length(dayopt)
                    % We need a group of units that is likely to be a pair across at least two days
                    if did==1
                        pairidx = recsesGood(Pairs(:,1)) == dayopt(did) & recsesGood(Pairs(:,2))==dayopt(did+1);
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
                        pairidx = recsesGood(Pairs(:,2)) == dayopt(did) & recsesGood(Pairs(:,1))==dayopt(did-1);
                        PairsTmp = Pairs(pairidx,:);
                        % Only use every 'unit' once --> take the highest scoring matches
                        [~,id1,~]=unique(PairsTmp(:,1),'stable');
                        PairsTmp = PairsTmp(id1,:);
                        [~,id1,~]=unique(PairsTmp(:,2),'stable');
                        PairsTmp = PairsTmp(id1,:);
                        Unit2TakeIdx = [Unit2TakeIdx; PairsTmp(:,2)];
                    end
                    Unit2TakeIdxAll = find(recsesGood == dayopt(did));
    
                    % Extract the part of the correlation matrix with these
                    % pairs
                    sortIdx = cell2mat(arrayfun(@(x) find(Unit2Take(Unit2TakeIdxAll) == x), Unit2Take(Unit2TakeIdx), 'uni', 0));
                    SessionCorrelation_Pair{did} = SessionCorrelations_perDay{dayopt(did)}(sortIdx,:);

                    % Remove =1 for the same unit (per def. 1)
                    for id = 1:length(Unit2TakeIdx)
                        SessionCorrelation_Pair{did}(id,ismember(Unit2TakeIdxAll,Unit2TakeIdx(id))) = nan;
                    end
    
                    % Normalize correlations to compare across recordings
                    SessionCorrelation_Pair{did} = (SessionCorrelation_Pair{did}-nanmedian(SessionCorrelation_Pair{did}(:)))./ ...
                        (quantile(SessionCorrelation_Pair{did}(:),0.95)-quantile(SessionCorrelation_Pair{did}(:),0.05));
                end

                % Get the correlation of the fingerprints
                FingerprintR = corr(SessionCorrelation_Pair{1},SessionCorrelation_Pair{2},'rows','pairwise');

                AllSessionCorrelations{did1,did2} = cat(2,SessionCorrelation_Pair{:})';
            end

            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1)-1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1)-1;

            % Save Fingerprint correlations
            FingerprintRAll(clusIdxD1All,clusIdxD2All) = FingerprintR;

            % Find rank
            FingerprintR(isnan(FingerprintR)) = 0; % should not participate to the rank
            [~,idx] = sort(FingerprintR,2,'descend');
            for c = 1:numel(clusIdxD1All)
                RankScoreAll(clusIdxD1All(c),clusIdxD2All(idx(c,:))) = 1:numel(clusIdxD2All);
            end

            % Find SigMask
            SigMask(clusIdxD1All,clusIdxD2All) = FingerprintR > quantile(FingerprintR,0.99,2);
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
                    imagesc(SessionCorrelations_fold1{did1}')
                    colormap(flipud(gray))
                    xlabel('Candidate Units to be matched')
                    ylabel(['Within day ' num2str(did1)])
                    title('First half of recording')
                    makepretty

                    subplot(nrows,2,rowcount*2)
                    imagesc(SessionCorrelations_fold2{did1}')
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
                        imagesc(AllSessionCorrelations{did1,did2}(clusIdxD{did},:)')
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

        %% Plot the rank

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
    end
end