function DataSizeParam = CalculateDuration(Dir2Use)
% Dir2Use = 'H:\Ongoing\EB019\'
tmp = dir(fullfile(Dir2Use,'UnitMatch','UnitMatch.mat'));

tmp = load(fullfile(tmp.folder,tmp.name));
TotalDur = 0;
TotalBytes = 0;
for id = 1:length(tmp.UMparam.AllRawPaths)

    meta = ReadMeta2(tmp.UMparam.AllRawPaths(id).folder);
    TotalDur = TotalDur+str2num(meta.fileTimeSecs);
    TotalBytes = TotalBytes + str2num(meta.fileSizeBytes);
end
nClus = sqrt(height(tmp.MatchTable));

DataSizeParam = struct;
DataSizeParam.TotalDur = TotalDur;
DataSizeParam.TotalBytes = TotalBytes;
DataSizeParam.nClus = nClus;
disp(['Total duration of recordings: ' num2str(round(TotalDur./60./60*10)/10) ' hours'])
disp(['Total size of recordings: ' num2str(round(TotalBytes./10^9*10)/10) ' GB'])
disp(['We extracted ' num2str(sqrt(height(tmp.MatchTable))) ' good clusters from this analysis'])