function [recordingNum, unitNum_0idx, rawWaveform, rawWaveformPath, bombcellPath, unitNum_bombcellRow] = um_getPathsAndInfo_forUID(thisUID, UniqueIDConversion)


    UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueID, 2), 1);
    UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;
    
    recordingNums = UniqueIDConversion.recsesAll(thisMatchTableIdx);
    
    for iRecordingNums = 1:size(recordingNums)
        thisRecording = recordingNums();
        unitNum_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
                UniqueIDConversion.UniqueID' == thisUID & ...
                UniqueIDConversion.recsesAll == thisRecording);
        unitNum_bombcellRow
        rawWaveformPath
    end


end