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
if any(ismember(MatchTable.Properties.VariableNames,'MatchProbUM')) % Indicates we compare UM to DNN
    % Ensure the same amount of GoodID
    % DUM:
    [MatchTableDUM, UniqueIDConversionDUM, UMparam] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID, 'MatchProb');
    % Original UM:
    [MatchTable, UniqueIDConversion, UMparam] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID, 'MatchProbUM');

    % Add them together
    UniqueIDConversion.UniqueIDConservativeDUM = UniqueIDConversionDUM.UniqueIDConservative;
    UniqueIDConversion.UniqueIDLiberalDUM = UniqueIDConversionDUM.UniqueIDLiberal;
    UniqueIDConversion.UniqueIDDUM = UniqueIDConversionDUM.UniqueID;
    if numel(UniqueIDConversion.UniqueID) ~= numel(UniqueIDConversion.UniqueIDDUM) 
        keyboard
    end

    % same for matchtable
    MatchTable.UID1ConservativeDUM = MatchTableDUM.UID1ConservativeDUM;
    MatchTable.UID2ConservativeDUM = MatchTableDUM.UID2ConservativeDUM;
    MatchTable.UID1LiberalDUM = MatchTableDUM.UID1LiberalDUM;
    MatchTable.UID2LiberalDUM = MatchTableDUM.UID2LiberalDUM;
    MatchTable.UID1DUM = MatchTableDUM.UID1DUM;
    MatchTable.UID2DUM = MatchTableDUM.UID2DUM;
    MatchTable.MatchProbDUM = MatchTableDUM.MatchProb; % Just to be certain, store this as the DUM as well 


else
    [MatchTable, UniqueIDConversion, UMparam] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam, StartUID);
end

%% Overwrite
disp('Saving table')
save(fullfile(tmpfile.folder,tmpfile.name),'MatchTable','UniqueIDConversion','-append')

return