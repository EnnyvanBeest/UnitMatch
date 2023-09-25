function DataSizeParam = CalculateDuration(Dir2Use, metaLocation)
% Dir2Use = 'H:\Ongoing\EB019\'
tmp = dir(fullfile(Dir2Use,'UnitMatch.mat'));

tmp = load(fullfile(tmp.folder,tmp.name));
TotalDur = 0;
TotalBytes = 0;
if nargin>1 && ~isempty(metaLocation)
    meta_path = metaLocation;
else
    meta_path = tmp.UMparam.AllRawPaths;
end
for id = 1:length(meta_path)
    if isstruct(meta_path{id})
        meta = ReadMeta2(meta_path{id}.folder);
    else
        meta_folder = dir(meta_path{id});
        meta = ReadMeta2(meta_folder.folder);
    end
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