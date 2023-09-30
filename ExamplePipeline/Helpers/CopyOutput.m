%% Path information
% Where to find and store data:
SaveDir = 'H:\Ongoing\'%'H:\SfN_2022'; %%'E:\Data\ResultsOngoing' %
CopyDir = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal'
% Spatial project %MiceOpt = {'EB001','EB002','EB003','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019','CB007','CB008','CB020','FT039','AV009','AV015','AV019','AV021','EB020','EB022','EB027','EB029','EB032','AV049'};%{'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019'};%{'EB019'};%% {'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019','CB007','CB008','CB020','FT039','AV009','AV015','AV019','AV021','EB020','EB022','EB027'}; % Add all mice you want to analyse
MiceOpt = {'AV015','EB014','EB019','AV008','AV021','FT033','FT039','AV009'};%{'AV009'};%;%{'CB017','CB018','FT033','AL032','AV008','CB016','JF067','JF078','JFAL35','EB014','EB019','CB007','CB008','CB020','FT039','AV009','AV015','AV019','AV021'};%{'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019'};%{'EB019'};%% {'EB001','EB002','EB003','EB004','EB005','EB006','EB007','EB008','EB009','EB010','EB011','EB012','EB013','EB014','EB017','EB016','EB019','CB007','CB008','CB020','FT039','AV009','AV015','AV019','AV021','EB020','EB022','EB027'}; % Chronic only mice
FromDate = datetime("2023-09-29 21:00:00");
for midx = 1:length(MiceOpt)
    fprintf('Reference %s...\n', MiceOpt{midx})
    % Identify all UM tables
    tmpfile = dir(fullfile(SaveDir, MiceOpt{midx},'*','*','UnitMatch', 'UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end
    for id = 1:length(tmpfile)
        if datetime(tmpfile(id).date) >FromDate
            FolderParts = strsplit(tmpfile(id).folder,filesep);
            idx = find(ismember(FolderParts,MiceOpt{midx}));
            if ~exist(fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2}))
                copyfile(fullfile(SaveDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2},'UnitMatch','*'),fullfile(CopyDir,FolderParts{idx},FolderParts{idx+1},FolderParts{idx+2}))
            end
        end
    end
end
