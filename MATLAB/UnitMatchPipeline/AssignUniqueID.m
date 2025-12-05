function [UniqueIDConversion, MatchTable, UMparam] = AssignUniqueID(SaveDir, StartUID)

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
if any(ismember(MatchTable.Properties.VariableNames,'NBProb18mice')) % Indicates we compare UM to DNN
    % Ensure the same amount of GoodID
    % DUM:
    display('DUM')
    [MatchTableDUM, UniqueIDConversionDUM, UMparam] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID, 'NBProb18mice');
    % Original UM:
    display('UM OG')

    [MatchTable, UniqueIDConversion, UMparam] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID, 'MatchProbNew');

    % Add them together
    UniqueIDConversion.UniqueIDConservativeDUM = UniqueIDConversionDUM.UniqueIDConservative;
    UniqueIDConversion.UniqueIDLiberalDUM = UniqueIDConversionDUM.UniqueIDLiberal;
    UniqueIDConversion.UniqueIDDUM = UniqueIDConversionDUM.UniqueID;
    if numel(UniqueIDConversion.UniqueID) ~= numel(UniqueIDConversion.UniqueIDDUM) 
        keyboard
    end

    % same for matchtable
    MatchTable.UID1ConservativeDUM = MatchTableDUM.UID1Conservative;
    MatchTable.UID2ConservativeDUM = MatchTableDUM.UID2Conservative;
    MatchTable.UID1LiberalDUM = MatchTableDUM.UID1Liberal;
    MatchTable.UID2LiberalDUM = MatchTableDUM.UID2Liberal;
    MatchTable.UID1DUM = MatchTableDUM.UID1;
    MatchTable.UID2DUM = MatchTableDUM.UID2;
    MatchTable.MatchProbDUM = MatchTableDUM.MatchProb; % Just to be certain, store this as the DUM as well 


else
    [MatchTable, UniqueIDConversion, UMparam] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID);
end

%% Overwrite
disp('Saving table')
save(fullfile(tmpfile.folder,tmpfile.name),'MatchTable','UniqueIDConversion','-append')

return