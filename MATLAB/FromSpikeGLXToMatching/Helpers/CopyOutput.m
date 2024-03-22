%% Path information
% Where to find and store data:
if 0 
SaveDir = 'H:\MatchingUnits\Output';%'H:\UnitMatch\'%'H:\SfN_2022'; %%'E:\Data\ResultsOngoing' %
CopyDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDaysNew\Stitched\'
% Spatial project %MiceOpt = {'EB001','EB002','EB003','EB005','EB006','EB007',EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019','CB007','CB008','CB020','FT039','AV009','AV015','AV019','AV021','EB020','EB022','EB027','EB029','EB032','AV049'};%{'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019'};%{'EB019'};%% {'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019','CB007','CB008','CB020','FT039','AV009','AV015','AV019','AV021','EB020','EB022','EB027'}; % Add all mice you want to analyse
% MiceOpt = {'CB016'};%{'JF067','FT033','AV008','JF084','JF078','JFAL35','AV009','EB014','EB019','FT039','AV015','AV021','AV049','AV019','EB036','EB037'};% %Chronic only
MiceOpt = {'AL032','AV008','CB016','EB019','JF067'}; %'AL032', Add all mice you want to analyze
% FromDate = datetime("2024-03-08 09:00:00");

FromDate = datetime("2024-03-16 14:00:00");
for midx = 1:length(MiceOpt)
    fprintf('Reference %s...\n', MiceOpt{midx})
    % Identify all UM tables
    tmpfile = dir(fullfile(SaveDir, MiceOpt{midx},'*','*','UnitMatch', 'UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end
    for id = 1:length(tmpfile)
        if datetime(tmpfile(id).date) > FromDate
%             AssignUniqueID_POSTUM(fullfile(tmpfile(id).folder,tmpfile(id).name));

            FolderParts = strsplit(tmpfile(id).folder,filesep);
            idx = find(ismember(FolderParts,MiceOpt{midx}));
            curFile = dir(fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2}));
            curFile(~[curFile.isdir]) = [];
            curFile(ismember({curFile(:).name},{'.','..'})) = [];

            if isempty(curFile) || curFile.date<FromDate
%                 ,'UnitMatch' change UMParam.SaveDir
                copyfile(fullfile(SaveDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2},'UnitMatch','*'),fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2},FolderParts{idx+3}))
                % now replace UMparam.SaveDIr
                load(fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2},FolderParts{idx+3},'UnitMatch.mat'),'UMparam')
                UMparam.SaveDir = fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2},FolderParts{idx+3});
                if ~isfield(UMparam,'RawDataPaths')
                    UMparam.RawDataPaths = UMparam.AllRawPaths;
                    rmfield(UMparam,'AllRawPaths')
                end
                save(fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2},FolderParts{idx+3},'UnitMatch.mat'),'UMparam','-append')
            end
% 
%             % Files to remove?
%             curFile = dir(fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2}));
%             curFile([curFile.isdir]) = [];
%             if ~isempty(curFile)
%                 arrayfun(@(X) delete(fullfile(curFile(X).folder,curFile(X).name)),1:length(curFile),'Uni',0)
%             end
     
        end
    end
end
end
%% Replace raw dir with server dirs
DataDir = {'\\znas\Subjects','\\zinu.cortexlab.net\Subjects','\\zaru.cortexlab.net\Subjects','\\zortex.cortexlab.net\Subjects'}%'\\znas\Subjects' %' Check DataDir2Use '\\128.40.198.18\Subjects',
MiceOpt = {'AL032','AV008','CB016','EB019','JF067'}; %'AL032', Add all mice you want to analyze

% Directory settings per mouse:
DataDir2Use = repmat(4,[1,length(MiceOpt)]);
DataDir2Use(ismember(MiceOpt,{'EB001','EB002','EB003','EB004','EB005','CB007','CB008','AL032'}))=1;
DataDir2Use(ismember(MiceOpt,{'JF067','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB015','AV008','AV009','CB020','FT039','EB017','EB018','CB016','FT033','CB017','CB018','CB020'}))=2;
DataDir2Use(ismember(MiceOpt,{'EB019','AV015','AV019','AV021','EB020','EB022','EB027','EB029','EB032','AV049','EB036','EB037'}))=3;

Dir2Check  = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDaysNew\NonStitched'

for midx = 1:length(MiceOpt)
    fprintf('Reference %s...\n', MiceOpt{midx})
    % Identify all UM tables
    tmpfile = dir(fullfile(Dir2Check, MiceOpt{midx},'*','*','UnitMatch', 'UnitMatch.mat'));
    if ~isempty(tmpfile)
        load(fullfile(tmpfile.folder,tmpfile.name),'UMparam');
        RawDatPath = UMparam.RawDataPaths;
        for rid = 1:length(RawDatPath)
            % if ~isdir(RawDatPath{rid})
            %     RawDatPath{rid} = dir(RawDatPath{rid});
            % end
            ServerPath = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},'**',RawDatPath{rid}.name));
            if length(ServerPath)==1
                UMparam.RawDataPaths{rid} = ServerPath;
                
            else
                keyboard
            end
        end
        % Make sure correct savedir
        Tmp = strsplit(UMparam.SaveDir,MiceOpt{midx});
        UMparam.SaveDir = fullfile(Dir2Check,MiceOpt{midx},Tmp{2});

        UMparam.AllRawPaths =  UMparam.RawDataPaths; % Duplicate
        save(fullfile(tmpfile.folder,tmpfile.name),'UMparam','-append');

        % ComputeFunctionalScores(UMparam.SaveDir,1,1)

    end



end
