function DataSetInfo = DataSetInfo(UMFiles,groupVector)
RecSes = 0;
if nargin<2 || ~exist('groupVector','var') || isempty(groupVector)
    groupVector = 1:length(UMFiles);
end
groups = unique(groupVector);
groupColor = distinguishable_colors(max(groups)+1);

days = cell(1, length(UMFiles));
deltaDays = cell(1, length(UMFiles));


% WIthin session check
EPosAndNeg = nan(2,length(UMFiles));
nGoodUnitsPerRec = cell(1,length(UMFiles));
nTotalUnitsPerRec = cell(1,length(UMFiles));
for midx = 1:length(UMFiles)

    %% Load data

    fprintf('Reference %s...\n', UMFiles{midx})

    tmpfile = dir(UMFiles{midx});
    if isempty(tmpfile)
        continue
    end

    %% Get some baseline information:
    fprintf('Loading the data...\n')
    tic
    load(fullfile(tmpfile.folder, tmpfile.name), 'UMparam', 'UniqueIDConversion');
    toc

    %% Sort out day situation
    if ~isstruct(UMparam.RawDataPaths{1})
        UMparam.RawDataPaths = cellfun(@(x) dir(x), UMparam.RawDataPaths, 'uni', 0);
    end
    days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    days{midx} = cell2mat(days{midx}) - days{midx}{1};
    deltaDays{midx} = days{midx} - days{midx}';
    
    %% Check if enough units?
    nclusGoodPerSess = arrayfun(@(X) sum(UniqueIDConversion.GoodID(UniqueIDConversion.recsesAll==X)),unique(UniqueIDConversion.recsesAll));
        nclusTotalPerSess = arrayfun(@(X) sum(UniqueIDConversion.recsesAll==X),unique(UniqueIDConversion.recsesAll));


    if all(nclusGoodPerSess<UMparam.minGoodUnits)
        continue
    end

    % Remove data for days that are having too little neurons
    if any(nclusGoodPerSess<UMparam.minGoodUnits)
        if median(nclusGoodPerSess)<UMparam.minGoodUnits
            continue
        end
    end

    RecSes = RecSes + numel(nclusGoodPerSess);
    nGoodUnitsPerRec{midx} = nclusGoodPerSess;
    nTotalUnitsPerRec{midx} = nclusTotalPerSess;
   
end



%% Info on dataset
DataSetInfo.RecSes = RecSes;
DataSetInfo.nGoodUnits = nGoodUnitsPerRec;
DataSetInfo.nTotalUnits = nTotalUnitsPerRec;

