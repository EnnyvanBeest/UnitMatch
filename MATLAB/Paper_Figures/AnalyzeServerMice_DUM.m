% 
SaveDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap'; %H:\Ongoing\'%'H:\SfN_2022'; %%'E:\Data\ResultsOngoing' %
% SaveDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays_KSChanMap\Stitched\'
% SaveDir = 'H:\UnitMatch\'
FromDate = datetime("2024-03-08 09:00:00");
AssignUnitDate = datetime("2026-01-17 10:21:00");

DUMFiles = cell(1,0); % Define your UMfolders here or use below:
UMFiles = cell(1,0);
groupvec = nan(1,0);
SpGLXV = cell(1,0);
if ~exist('UMFiles') || isempty(DUMFiles) % When using the example pipeline this may be useful:
    MiceOpt = dir(SaveDir);
    MiceOpt = arrayfun(@(X) X.name,MiceOpt,'Uni',0);
    MiceOpt(ismember(MiceOpt,{'.','..'})) = [];

    for midx = 1:numel(MiceOpt)
        fprintf('Reference %s...\n', MiceOpt{midx})
        % Identify all UM tables
        % tmpfile = dir(fullfile(SaveDir, MiceOpt{midx}, 'UnitMatch.mat'));

        tmpfile = dir(fullfile(SaveDir, MiceOpt{midx},'*','*','DeepUnitMatch', 'UnitMatch.mat'));

        if isempty(tmpfile)
            continue
        end
        for id = 1:length(tmpfile)
            % if datetime(tmpfile(id).date) > FromDate % && any(cell2mat(cellfun(@(X) any(strfind(fullfile(tmpfile(id).folder,tmpfile(id).name),X)),UMFiles2Take,'Uni',0)))
                
            if datetime(tmpfile(id).date) < AssignUnitDate
                AssignUniqueID(fullfile(tmpfile(id).folder)) % REDO
            end
            % Check that these data are not too noisy

                % load(fullfile(tmpfile(id).folder,tmpfile(id).name),'UMparam')
                % for rid = 1:numel(UMparam.RawDataPaths)
                %    meta = ReadMeta2(UMparam.RawDataPaths{rid}.folder);
                %    SpGLXV = {SpGLXV{:} meta.appVersion};
                % end

                % Check some stuff across DUM folder and UM folder
                


                %             FolderParts = strsplit(tmpfile(id).folder,filesep);
                %             idx = find(ismember(FolderParts,MiceOpt{midx}));
                DUMFiles = cat(2,DUMFiles,fullfile(tmpfile(id).folder,tmpfile(id).name));

                % also store the UM files
                UMFiles = cat(2,UMFiles,strrep(fullfile(tmpfile(id).folder,tmpfile(id).name),'DeepUnitMatch','UnitMatch'));

                groupvec = cat(2,groupvec,midx);
            % else
                % keyboard
            % end


        end
    end
    close all
end

Info  = DataSetInfo(DUMFiles)
Info.RecSes
nanmean(cat(1,Info.nGoodUnits{:})./cat(1,Info.nTotalUnits{:}).*100)
nanstd(cat(1,Info.nGoodUnits{:})./cat(1,Info.nTotalUnits{:}).*100)



InfoUM  = DataSetInfo(UMFiles)
InfoUM.RecSes
nanmean(cat(1,InfoUM.nGoodUnits{:})./cat(1,InfoUM.nTotalUnits{:}).*100)
nanstd(cat(1,InfoUM.nGoodUnits{:})./cat(1,InfoUM.nTotalUnits{:}).*100)

%%
AUCScoresDUM = AUCExtract(DUMFiles)
AUCScoresUM = AUCExtract(UMFiles)


%%
[unitPresence, unitProbaMatch, days, EPosAndNeg, DataSetInfo, pSinceAppearance, popCorr_Uni, popAUC_Uni] = summaryMatchingPlots(DUMFiles,{'UID1','UID1DUM'},groupvec,1,UMFiles)

%%
UIDTypes = {'UM','DUM'};

% UIDTypes = {'UM','DUM','OriUM'}; % When comparing to original UM paper
% (in case anyone asks)
UIDCols = [0 0.7 0.2; 0 0 0; 0.7 0.2 0];
ISI = DataSetInfo.ISI(:,:,1:numel(UIDTypes));
pSinceAppearance = pSinceAppearance(:,:,1:numel(UIDTypes));

%
avgAcrossMice = arrayfun(@(X) squeeze(nanmean(pSinceAppearance(:,groupvec==X,:),2)),unique(groupvec),'Uni',0);
avgAcrossMice = cat(3,avgAcrossMice{:});

avgISIAcrossMice = arrayfun(@(X) squeeze(nanmean(ISI(:,groupvec==X,:),2)),unique(groupvec),'Uni',0);
avgISIAcrossMice = cat(3,avgISIAcrossMice{:});
% Mixed effects model.
% Test UIDType, DaysBin (categorical), and interaction with mouse as random intercept.
nBins = size(avgAcrossMice, 1);
nMice = size(avgAcrossMice, 3);
nUid = numel(UIDTypes);
[binIdx, dsIdx, uidIdx] = ndgrid(1:nBins, 1:nMice, 1:nUid);
resp = avgAcrossMice(:);
valid = ~isnan(resp);
tbl = table;
epsVal = 1e-3;
respT = log((resp + epsVal) ./ (1 - resp + epsVal));
tbl.Response = respT(valid);
tbl.ISI = avgISIAcrossMice(valid);
tbl.DaysBin = categorical(binIdx(valid));
tbl.UIDType = categorical(uidIdx(valid), 1:nUid, UIDTypes);

%% Statistical model:
dayGroups = findgroups(tbl.DaysBin);
uidGroups = findgroups(tbl.UIDType);
comboCounts = accumarray([dayGroups, uidGroups], 1);
keepDays = all(comboCounts > 0, 2);
daysAll = categories(tbl.DaysBin);
keepDayCats = daysAll(keepDays);
tbl = tbl(ismember(tbl.DaysBin, keepDayCats), :);
% tbl.DaysBin = removecats(tbl.DaysBin);
% tbl.UIDType = removecats(tbl.UIDType);
lm = fitlm(tbl, 'Response ~ UIDType*DaysBin');
anova(lm)
effectNames = lm.CoefficientNames;
daysCats = categories(tbl.DaysBin);
dayCatNums = str2double(daysCats);
uidCats = categories(tbl.UIDType);
nDays = numel(daysCats);
nUid = numel(uidCats);
pVals_UID = nan(nUid, nUid, size(avgAcrossMice,1));
diffEst_UID = nan(nUid, nUid, size(avgAcrossMice,1));
pHolm_UID = nan(nUid, nUid, size(avgAcrossMice,1));
ttestvals = nan(nUid,nUid,size(avgISIAcrossMice,1));

for dci = 1:nDays
    dayName = daysCats{dci};
    dayNum = dayCatNums(dci);
    if isnan(dayNum)
        continue
    end
    for uid1 = 1:nUid
        v1 = design_vec(uidCats{uid1}, dayName, effectNames, uidCats{1}, daysCats{1});
        for uid2 = 1:nUid
            if uid1 == uid2
                continue
            end
            v2 = design_vec(uidCats{uid2}, dayName, effectNames, uidCats{1}, daysCats{1});
            contrast = v1 - v2;
            [p, ~, ~] = coefTest(lm, contrast);
            pVals_UID(uid1, uid2, dayNum) = p;
            diffEst_UID(uid1, uid2, dayNum) = contrast * lm.Coefficients.Estimate;
        end
    end
    idxTri = find(triu(ones(nUid), 1));
    pVec = pVals_UID(:, :, dayNum);
    pAdj = holm_adjust(pVec(idxTri));
    pHolm = nan(nUid);
    pHolm(idxTri) = pAdj;
    pHolm_UID(:, :, dayNum) = pHolm;
    ttestvals(:,:,dayNum) = signrank(squeeze(avgAcrossMice(dayNum,1,:)),squeeze(avgAcrossMice(dayNum,2,:)))

end
posthoc = struct();
posthoc.UIDTypes = uidCats;
posthoc.DaysBin = daysCats;
posthoc.pVals = pVals_UID;
posthoc.pHolm = pHolm_UID;
posthoc.diffEst = diffEst_UID;
posthoc.ttestvals = ttestvals;

disp(posthoc)

%
clear h
figure
subplot(1,2,1)
title('Proportion track')
hold on
for uid1 = 1:nUid
    nonnanNr = sum(~isnan(avgAcrossMice(:,uid1,:)),3)+1; %nr data +1
    h(uid1) = errorbar(1:size(avgAcrossMice,1),nanmean(avgAcrossMice(:,uid1,:),3),nanstd(avgAcrossMice(:,uid1,:),[],3)./sqrt(nonnanNr-1),'linestyle','-','color',UIDCols(uid1,:));

    % shadedErrorBar(1:size(pSinceAppearance,1),nanmean(pSinceAppearance,2),nanstd(pSinceAppearance,[],2)./sqrt(nonnanNr-1),'transparent',1)


    for uid2 = 1:nUid
        if uid1 >= uid2
            continue
        end
        if any(squeeze(ttestvals(uid1,uid2,:)<0.05))
        plot(find(squeeze(ttestvals(uid1,uid2,:)<0.05)),...
            repmat(max(get(gca,'ylim'))+0.1,sum(squeeze(ttestvals(uid1,uid2,:)<0.05),1)),'*','color',nanmean(UIDCols([uid1,uid2],:),1));
       
        end
        ylim([0 max(get(gca,'ylim'))+0.05])
    end
end
xlabel('delta Days (>)')
ylabel('P(track)')
makepretty
offsetAxes
ylim([0 0.6])
legend([h(:)],UIDTypes)


%% Same for ISI
tbl = table;
epsVal = 1e-3;
respT = log((resp + epsVal) ./ (1 - resp + epsVal));
tbl.Response = respT(valid);
tbl.ISI = avgISIAcrossMice(valid);
tbl.DaysBin = categorical(binIdx(valid));
tbl.UIDType = categorical(uidIdx(valid), 1:nUid, UIDTypes);
tblISI = tbl(~isnan(tbl.ISI), :);
% tblISI.DaysBin = removecats(tblISI.DaysBin);
% tblISI.UIDType = removecats(tblISI.UIDType);
lmISI = fitlm(tblISI, 'ISI ~ UIDType*DaysBin');
anova(lmISI)
effectNames = lmISI.CoefficientNames;
daysCats = categories(tblISI.DaysBin);
dayCatNums = str2double(daysCats);
uidCats = categories(tblISI.UIDType);
nDays = numel(daysCats);
nUid = numel(uidCats);
pVals_UID = nan(nUid, nUid, size(avgISIAcrossMice,1));
diffEst_UID = nan(nUid, nUid, size(avgISIAcrossMice,1));
pHolm_UID = nan(nUid, nUid, size(avgISIAcrossMice,1));

ttestvals = nan(nUid,nUid,size(avgISIAcrossMice,1));
for dci = 1:nDays
    dayName = daysCats{dci};
    dayNum = dayCatNums(dci);
    if isnan(dayNum)
        continue
    end
    for uid1 = 1:nUid
        v1 = design_vec(uidCats{uid1}, dayName, effectNames, uidCats{1}, daysCats{1});
        for uid2 = 1:nUid
            if uid1 == uid2
                continue
            end
            v2 = design_vec(uidCats{uid2}, dayName, effectNames, uidCats{1}, daysCats{1});
            contrast = v1 - v2;
            [p, ~, ~] = coefTest(lmISI, contrast);
            pVals_UID(uid1, uid2, dayNum) = p;
            diffEst_UID(uid1, uid2, dayNum) = contrast * lmISI.Coefficients.Estimate;
        end
    end
    idxTri = find(triu(ones(nUid), 1));
    pVec = pVals_UID(:, :, dayNum);
    pAdj = holm_adjust(pVec(idxTri));
    pHolm = nan(nUid);
    pHolm(idxTri) = pAdj;
    pHolm_UID(:, :, dayNum) = pHolm;
    ttestvals(:,:,dayNum) = signrank(squeeze(avgISIAcrossMice(dayNum,1,:)),squeeze(avgISIAcrossMice(dayNum,2,:)))
end
posthoc = struct();
posthoc.UIDTypes = uidCats;
posthoc.DaysBin = daysCats;
posthoc.pVals = pVals_UID;
posthoc.pHolm = pHolm_UID;
posthoc.diffEst = diffEst_UID;
posthoc.ttestvals = ttestvals;

disp(posthoc)



%
subplot(1,2,2)
title('ISI AUC')
hold on
for uid1 = 1:nUid
    nonnanNr = sum(~isnan(ISI(:,:,uid1)),2)+1; %nr data +1
    h(uid1) = errorbar(1:size(avgISIAcrossMice,1),nanmean(avgISIAcrossMice(:,uid1,:),3),nanstd(avgISIAcrossMice(:,uid1,:),[],3)./sqrt(nonnanNr-1),'linestyle','-','color',UIDCols(uid1,:));
  
    % shadedErrorBar(1:size(pSinceAppearance,1),nanmean(pSinceAppearance,2),nanstd(pSinceAppearance,[],2)./sqrt(nonnanNr-1),'transparent',1)
   

    for uid2 = 1:nUid
        if uid1 >= uid2
            continue
        end
         if any(squeeze(ttestvals(uid1,uid2,:)<0.05))
        plot(find(squeeze(ttestvals(uid1,uid2,:)<0.05)),...
            repmat(max(get(gca,'ylim'))+0.1,sum(squeeze(ttestvals(uid1,uid2,:)<0.05),1)),'*','color',nanmean(UIDCols([uid1,uid2],:),1));
        end
        ylim([0.5 max(get(gca,'ylim'))+0.05])
    end
end
xlabel('delta Days (>)')
ylabel('ISI AUC')
makepretty
offsetAxes
legend([h(:)],UIDTypes)
% 

%%
summaryFunctionalPlots_Part2(DUMFiles, groupvec, 0)

function pAdj = holm_adjust(pVals)
pAdj = nan(size(pVals));
pFlat = pVals(:);
valid = ~isnan(pFlat);
if ~any(valid)
    return
end
idxValid = find(valid);
[pSorted, order] = sort(pFlat(idxValid));
m = numel(pSorted);
adjSorted = cummax((m - (0:m-1)') .* pSorted);
adjSorted = min(adjSorted, 1);
pAdjValid = zeros(m, 1);
pAdjValid(order) = adjSorted;
pFlat(idxValid) = pAdjValid;
pAdj(:) = pFlat;
end

function v = design_vec(uidName, dayName, effectNames, uidRef, dayRef)
v = zeros(1, numel(effectNames));
idxIntercept = strcmp(effectNames, '(Intercept)');
v(idxIntercept) = 1;
if ~strcmp(uidName, uidRef)
    idxUID = strcmp(effectNames, ['UIDType_' uidName]);
    v(idxUID) = 1;
end
if ~strcmp(dayName, dayRef)
    idxDay = strcmp(effectNames, ['DaysBin_' dayName]);
    v(idxDay) = 1;
end
if ~strcmp(uidName, uidRef) && ~strcmp(dayName, dayRef)
    term1 = ['UIDType_' uidName ':DaysBin_' dayName];
    term2 = ['DaysBin_' dayName ':UIDType_' uidName];
    idxInt = strcmp(effectNames, term1) | strcmp(effectNames, term2);
    v(idxInt) = 1;
end
end


