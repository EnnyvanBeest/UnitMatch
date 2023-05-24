% function ComputeFunctionalScores(SaveDir)

TmpFile = matfile(fullfile(SaveDir,'UnitMatch','UnitMatch.mat')); % Access saved file
UMparam = TmpFile.UMparam; % Extract parameters

MatchTable = TmpFile.MatchTable; % Load Matchtable

% Extract groups
WithinIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
MatchIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
NonMatchIdx = find((MatchTable.UID1 ~= MatchTable.UID2)); % Not the same unit

% Extract cluster information
UniqueIDConversion = TmpFile.UniqueIDConversion;
GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
OriID = UniqueIDConversion.OriginalClusID(GoodId);
recses = UniqueIDConversion.recsesAll(GoodId);
AllKSDir = UMparam.KSDir; %original KS Dir
nclus = length(UniqueID);


%% Cross-correlation ROC?
figure('name','Cross-correlation')
subplot(2,2,1)
bins = -1:0.1:1;
Vector = [-1+0.1/2:0.1:1-0.1/2];
hw = histcounts(MatchTable.FingerprintCor(WithinIdx),bins)./length(WithinIdx);
hm = histcounts(MatchTable.FingerprintCor(MatchIdx),bins)./length(MatchIdx);
hn = histcounts(MatchTable.FingerprintCor(NonMatchIdx),bins)./length(NonMatchIdx);
plot(Vector,hw,'color',[0.5 0.5 0.5])
hold on
plot(Vector,hm,'color',[0 0.5 0])
plot(Vector,hn,'color',[0 0 0])
xlabel('Cross-correlation Fingerprint')
ylabel('Proportion|Group')
legend('i=j; within recording','matches','non-matches','Location','best')
makepretty

subplot(2,2,2)
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(NonMatchIdx))];
scores = [MatchTable.FingerprintCor(MatchIdx)', MatchTable.FingerprintCor(NonMatchIdx)'];
[X,Y,~,AUC1] = perfcurve(labels,scores,1);
h(1) = plot(X,Y,'color',[0 0.25 0]);
hold all
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
scores = [MatchTable.FingerprintCor(MatchIdx)', MatchTable.FingerprintCor(WithinIdx)'];
[X,Y,~,AUC2] = perfcurve(labels,scores,1);
h(2) = plot(X,Y,'color',[0 0.5 0]);
labels = [ones(1,numel(WithinIdx)), zeros(1,numel(NonMatchIdx))];
scores = [MatchTable.FingerprintCor(WithinIdx)', MatchTable.FingerprintCor(NonMatchIdx)'];
[X,Y,~,AUC3] = perfcurve(labels,scores,1);
h(3) = plot(X,Y,'color',[0.25 0.25 0.25]);

plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend([h(:)],'Match vs No Match','Match vs Within','Within vs No Match','Location','best')
title(sprintf('AUC: %.3f, %.3f, %.3f', AUC1,AUC2,AUC3))
makepretty
drawnow %Something to look at while ACG calculations are ongoing

if ~any(ismember(MatchTable.Properties.VariableNames,'ACGCorr')) % If it already exists in table, skip this entire thing
    % Load SP
    disp('Loading spike information...')
    nKSFiles = length(AllKSDir);
    sp = cell(1,nKSFiles);
    for did = 1:nKSFiles
        tmp = matfile(fullfile(AllKSDir{did},'PreparedData.mat'));
        sptmp = tmp.sp;
        clear tmp

        % Only keep parameters used
        sp{did}.st = sptmp.st;
        sp{did}.spikeTemplates = sptmp.spikeTemplates;
        sp{did}.spikeAmps = sptmp.spikeAmps;
        % Replace recsesid with subsesid
        sp{did}.RecSes = repmat(did,size(sp{did}.st));
    end
    clear sptmp


    % Add all spikedata in one spikes struct - can be used for further analysis
    sp = [sp{:}];
    spnew = struct;
    fields = fieldnames(sp(1));
    for fieldid=1:length(fields)
        try
            eval(['spnew.' fields{fieldid} '= cat(1,sp(:).' fields{fieldid} ');'])
        catch ME
            if strcmp(ME.message,'Out of memory.')
                eval(['spnew.' fields{fieldid} ' = sp(1).' fields{fieldid} ';'])
                for tmpid = 2:length(sp)
                    eval(['spnew.' fields{fieldid} ' = cat(1,spnew.' fields{fieldid} ', sp(tmpid).' fields{fieldid} ');'])
                end
            else
                eval(['spnew.' fields{fieldid} '= cat(2,sp(:).' fields{fieldid} ');'])
            end
        end
    end
    sp = spnew;
    clear spnew


    %% Compute ACG and correlate them between units
    % This is very time consuming
    disp('Computing ACG, this will take some time...')
    tvec =  -UMparam.ACGduration/2:UMparam.ACGbinSize:UMparam.ACGduration/2;
    ACGMat = nan(length(tvec),2,nclus);
    parfor clusid = 1:nclus
        for cv = 1:2
            idx1=find(sp.spikeTemplates == OriID(clusid) & sp.RecSes == recses(clusid));
            if ~isempty(idx1)
                if cv==1
                    idx1 = idx1(1:floor(length(idx1)/2));
                else
                    idx1 = idx1(ceil(length(idx1)/2):end);
                end

                % compute ACG
                [ccg, t] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
                    ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
                ACGMat(:,cv,clusid) = ccg(:, 1, 1);
            end
        end
    end

    %% Correlation between ACG
    ACGCorr = corr(squeeze(ACGMat(:,1,:)),squeeze(ACGMat(:,2,:)));
    MatchTable.ACGCorr = ACGCorr(:);
    TmpFile.Properties.Writable = true;
    TmpFile.MatchTable = MatchTable; % Overwrite
    TmpFile.Properties.Writable = false;
end

%% Plot ACG
subplot(2,2,3)
bins = -1:0.1:1;
Vector = [-1+0.1/2:0.1:1-0.1/2];
hw = histcounts(MatchTable.ACGCorr(WithinIdx),bins)./length(WithinIdx);
hm = histcounts(MatchTable.ACGCorr(MatchIdx),bins)./length(MatchIdx);
hn = histcounts(MatchTable.ACGCorr(NonMatchIdx),bins)./length(NonMatchIdx);
plot(Vector,hw,'color',[0.5 0.5 0.5])
hold on
plot(Vector,hm,'color',[0 0.5 0])
plot(Vector,hn,'color',[0 0 0])
xlabel('Autocorrelogram Correlation')
ylabel('Proportion|Group')
legend('i=j; within recording','matches','non-matches','Location','best')
makepretty

subplot(2,2,4)
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(NonMatchIdx))];
scores = [MatchTable.ACGCorr(MatchIdx)', MatchTable.ACGCorr(NonMatchIdx)'];
[X,Y,~,AUC1] = perfcurve(labels,scores,1);
h(1) = plot(X,Y,'color',[0 0.25 0]);
hold all
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
scores = [MatchTable.ACGCorr(MatchIdx)', MatchTable.ACGCorr(WithinIdx)'];
[X,Y,~,AUC2] = perfcurve(labels,scores,1);
h(2) = plot(X,Y,'color',[0 0.5 0]);
labels = [ones(1,numel(WithinIdx)), zeros(1,numel(NonMatchIdx))];
scores = [MatchTable.ACGCorr(WithinIdx)', MatchTable.ACGCorr(NonMatchIdx)'];
[X,Y,~,AUC3] = perfcurve(labels,scores,1);
h(3) = plot(X,Y,'color',[0.25 0.25 0.25]);

plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend([h(:)],'Match vs No Match','Match vs Within','Within vs No Match','Location','best')
title(sprintf('Autocorrelogram AUC: %.3f, %.3f, %.3f', AUC1,AUC2,AUC3))
makepretty