TotalDur = 0
TotalBytes = 0;
for id = 1:length(tmp.UMparam.AllRawPaths)

    meta = ReadMeta2(tmp.UMparam.AllRawPaths(id).folder);

    TotalDur = TotalDur+str2num(meta.fileTimeSecs)
TotalBytes = TotalBytes + str2num(meta.fileSizeBytes);
end