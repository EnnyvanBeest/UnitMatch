function [UniqueIDConversion, MatchTable] = AssignUniqueID(SaveDir, StartUID)

if nargin<2
    StartUID = 1;
end
disp('Loading match table to assign Unique IDs...')
tmpfile = dir(fullfile(SaveDir,'UnitMatch.mat'));
load(fullfile(tmpfile.folder,tmpfile.name))
if exist('TmpFile', 'var')
    UniqueIDConversion = TmpFile.UniqueIDConversion;
    MatchTable = TmpFile.MatchTable;
    UMparam = TmpFile.UMparam;
end
if ~isfield(UMparam,'UseDatadrivenProbThrs')
    UMparam.UseDatadrivenProbThrs = 0;
end
%% Subfunction in case you want to use this on a subset of the data/table
[MatchTable, UniqueIDConversion] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID);

%% Overwrite
disp('Saving table')
save(fullfile(tmpfile.folder,tmpfile.name),'MatchTable','UniqueIDConversion','-append')

return