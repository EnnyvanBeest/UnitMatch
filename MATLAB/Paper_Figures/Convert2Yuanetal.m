function Convert2Yuanetal(PathNames,ToPath)
paramBC = bc_qualityParamValuesForUnitMatch;

for did = 1:length(PathNames)
    if ~isfolder(fullfile(ToPath,['D' num2str(did)]))
        mkdir(fullfile(ToPath,['D' num2str(did)]))
    end
    Filetmp = dir(fullfile(PathNames{did},'**','channel_map.npy'));
    copyfile(fullfile(Filetmp.folder,Filetmp.name),fullfile(ToPath,['D' num2str(did)],Filetmp.name));

    Filetmp = dir(fullfile(PathNames{did},'**','channel_positions.npy'));
    copyfile(fullfile(Filetmp.folder,Filetmp.name),fullfile(ToPath,['D' num2str(did)],Filetmp.name));

    Filetmp = dir(fullfile(PathNames{did},'**','cluster_KSLabel.tsv'));
    copyfile(fullfile(Filetmp(1).folder,Filetmp(1).name),fullfile(ToPath,['D' num2str(did)],Filetmp(1).name));
    clusinfo = tdfread(fullfile(Filetmp(1).folder,Filetmp(1).name));


    % Replace with bombcell good units to compare to UnitMatch
    d = dir(fullfile(PathNames{did}, '**', 'templates._bc_qMetrics.parquet'));
    qMetricsPath = d.folder;
    [~, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(qMetricsPath);
    unitType = bc_getQualityUnitType(paramBC, qMetric);

    Filetmp = dir(fullfile(PathNames{did},'**','RawWaveforms'));
    Filetmp(ismember({Filetmp(:).name},{'.','..'})) = [];
    AllSpks = nan(numel(unitType),384,82);
    for uid = 1:length(Filetmp)
        clusid = strsplit(Filetmp(uid).name,'Unit');
        clusid = strsplit(clusid{2},'_');
        clusid = find(clusinfo.cluster_id==str2num(clusid{1}));
        tmp = readNPY(fullfile(Filetmp(uid).folder,Filetmp(uid).name));
        tmp = nanmean(tmp,3)';
        AllSpks(clusid,:,:) = tmp;
        if any(isnan(tmp(:))) & unitType(uid)==1
            unitType(uid) = 0;
        end
    end
    save(fullfile(ToPath,['D' num2str(did)],'Bombcellgood.mat'),'unitType');

    writeNPY(AllSpks,fullfile(ToPath,['D' num2str(did)],'ksproc_mean_waveforms.npy'))
end
return