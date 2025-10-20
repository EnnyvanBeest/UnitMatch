% 
SaveDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap'; %H:\Ongoing\'%'H:\SfN_2022'; %%'E:\Data\ResultsOngoing' %
% SaveDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays_KSChanMap\Stitched\'
% SaveDir = 'H:\UnitMatch\'
FromDate = datetime("2024-03-08 09:00:00");
AssignUnitDate = datetime("2025-10-15 14:21:00");

UMFiles = cell(1,0); % Define your UMfolders here or use below:
groupvec = nan(1,0);
SpGLXV = cell(1,0);
if ~exist('UMFiles') || isempty(UMFiles) % When using the example pipeline this may be useful:
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


                %             FolderParts = strsplit(tmpfile(id).folder,filesep);
                %             idx = find(ismember(FolderParts,MiceOpt{midx}));
                UMFiles = cat(2,UMFiles,fullfile(tmpfile(id).folder,tmpfile(id).name));
                groupvec = cat(2,groupvec,midx);
            % else
                % keyboard
            % end


        end
    end
    close all
end

Info  = DataSetInfo(UMFiles)
Info.RecSes
nanmean(cat(1,Info.nGoodUnits{:})./cat(1,Info.nTotalUnits{:}).*100)
nanstd(cat(1,Info.nGoodUnits{:})./cat(1,Info.nTotalUnits{:}).*100)


AUCExtract(UMFiles)


summaryMatchingPlots(UMFiles,{'UID1','UID1DUM'},groupvec,1)

%
summaryFunctionalPlots_Part2(UMFiles, groupvec, 0)


