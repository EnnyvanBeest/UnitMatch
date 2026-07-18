function [MatchTable, UMparam, UniqueIDConversion, WaveformInfo] = loadPythonOutputForMatlab(SaveDir)
%LOADPYTHONOUTPUTFORMATLAB Load Python batch outputs directly in MATLAB.
%
% This loader reads the sidecar CSV and optional Python-side files generated
% by the batch script and constructs the minimal structures expected by the
% MATLAB functions ComputeFunctionalScores and PlotUnitsOnProbe when
% UnitMatch.mat is not available.

    if nargin < 1 || isempty(SaveDir)
        error('SaveDir must be provided.');
    end

    SaveDir = char(SaveDir);
    if ~exist(SaveDir, 'dir')
        error('SaveDir does not exist: %s', SaveDir);
    end

    csvPath = fullfile(SaveDir, 'MatchTable.csv');
    if exist(csvPath, 'file')
        T = readtable(csvPath);
    else
        T = table();
    end

    if istable(T) && ~isempty(T)
        renameMap = containers.Map({'RecSes 1','RecSes 2','UMProbabilities','UM Probabilities','UID orig 1','UID orig 2','UID Lib 1','UID Lib 2','UID 1','UID 2','UID Cons 1','UID Cons 2','refpop_correlation','refpop_correlations','ISI_correlations','FR_diff','natim_correlation','natim_correlations'}, ...
            {'RecSes1','RecSes2','MatchProb','MatchProb','UID1','UID2','UID1','UID2','UID1','UID2','UID1','UID2','refPopCorr','refPopCorr','ISICorr','FRDiff','natImRespCorr','natImRespCorr'});

        varNames = T.Properties.VariableNames;
        for idx = 1:numel(varNames)
            varName = char(varNames{idx});
            if isKey(renameMap, varName)
                newName = char(renameMap(varName));
                if ~strcmp(varName, newName)
                    if ~ismember(newName, varNames(~ismember(varNames, varName)))
                        T.Properties.VariableNames{idx} = newName;
                    end
                end
            end
        end

        if ~ismember('MatchProb', T.Properties.VariableNames)
            T.MatchProb = nan(height(T), 1);
        end
        if ~ismember('TotalScore', T.Properties.VariableNames)
            T.TotalScore = T.MatchProb;
        end
        if ~ismember('EucledianDistance', T.Properties.VariableNames)
            T.EucledianDistance = nan(height(T), 1);
        end
        if ~ismember('UID1', T.Properties.VariableNames)
            T.UID1 = nan(height(T), 1);
        end
        if ~ismember('UID2', T.Properties.VariableNames)
            T.UID2 = nan(height(T), 1);
        end
        if ~ismember('RecSes1', T.Properties.VariableNames)
            T.RecSes1 = nan(height(T), 1);
        end
        if ~ismember('RecSes2', T.Properties.VariableNames)
            T.RecSes2 = nan(height(T), 1);
        end
        if ~ismember('ID1', T.Properties.VariableNames)
            T.ID1 = nan(height(T), 1);
        end
        if ~ismember('ID2', T.Properties.VariableNames)
            T.ID2 = nan(height(T), 1);
        end

        matlabCols = {'ID1','ID2','RecSes1','RecSes2','UID1','UID2','UIDCons1','UIDCons2','UIDLib1','UIDLib2','MatchProb','TotalScore','EucledianDistance','refPopCorr','ISICorr','FRDiff','natImRespCorr'};
        for col = matlabCols
            colName = char(col);
            if ~ismember(colName, T.Properties.VariableNames)
                T.(colName) = nan(height(T), 1);
            end
        end
        T = T(:, matlabCols);
    else
        T = table(nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), ...
            'VariableNames', {'ID1','ID2','RecSes1','RecSes2','UID1','UID2','UIDCons1','UIDCons2','UIDLib1','UIDLib2','MatchProb','TotalScore','EucledianDistance','refPopCorr','ISICorr','FRDiff','natImRespCorr'});
    end

    MatchTable = T;

    umparamPath = fullfile(SaveDir, 'UMparam.pickle');
    clusInfoPath = fullfile(SaveDir, 'ClusInfo.pickle');

    UMparam = struct();
    UMparam.GoodUnitsOnly = true;
    UMparam.ProbabilityThreshold = 0.5;
    UMparam.binsz = 0.01;

    if exist(umparamPath, 'file')
        try
            umparamPy = readPythonPickle(umparamPath);
            if isstruct(umparamPy)
                if isfield(umparamPy, 'KSDir') && ~isempty(umparamPy.KSDir)
                    UMparam.KSDir = umparamPy.KSDir;
                elseif isfield(umparamPy, 'KS_dirs') && ~isempty(umparamPy.KS_dirs)
                    UMparam.KSDir = umparamPy.KS_dirs;
                elseif isfield(umparamPy, 'KSdirs') && ~isempty(umparamPy.KSdirs)
                    UMparam.KSDir = umparamPy.KSdirs;
                end
            end
        catch
        end
    end


    UniqueIDConversion.UniqueID = T.UID1(:);
    UniqueIDConversion.UniqueIDConservative = T.UIDCons1(:);
    UniqueIDConversion.UniqueIDLiberal = T.UIDLib1(:);
    UniqueIDConversion.OriginalClusID = T.ID1(:);
    UniqueIDConversion.recsesAll = T.RecSes1(:);
    UniqueIDConversion = struct2table(UniqueIDConversion);

    UniqueIDConversion = unique(UniqueIDConversion,'stable');

    
    goodID = true(1, numel(UniqueIDConversion.UniqueID));
    if exist(clusInfoPath, 'file')
        try
            clusInfo = readPythonPickle(clusInfoPath);
            if isstruct(clusInfo) && isfield(clusInfo, 'good_units') && ~isempty(clusInfo.good_units)
                goodUnits = clusInfo.good_units;
                goodID = false(1, numel(UniqueIDConversion.UniqueID));
                for idx = 1:numel(goodUnits)
                    if isempty(goodUnits{idx})
                        continue
                    end
                    if iscell(goodUnits{idx})
                        unitVals = cell2mat(goodUnits{idx});
                    else
                        unitVals = goodUnits{idx};
                    end
                    unitVals = unitVals(:)';
                    goodID(ismember(UniqueIDConversion.OriginalClusID, unitVals) | ismember(UniqueIDConversion.UniqueID, unitVals)) = true;
                end
                if ~any(goodID)
                    goodID = true(1, numel(UniqueIDConversion.UniqueID));
                end
            end
        catch
            goodID = true(1, numel(UniqueIDConversion.UniqueID));
        end
    end
    UniqueIDConversion.GoodID = goodID';

    WaveformInfo = struct();
    wfPath = fullfile(SaveDir, 'WaveformInfo.npz');
    if exist(wfPath, 'file')
        try
            WaveformInfo = readPythonNpz(wfPath);
        catch
            WaveformInfo = struct();
        end
    end
    if ~isfield(WaveformInfo, 'ProjectedLocation')
        WaveformInfo.ProjectedLocation = zeros(3, numel(UniqueIDConversion.UniqueID), 1);
    end
    if ~isfield(WaveformInfo, 'ProjectedWaveform')
        WaveformInfo.ProjectedWaveform = [];
    end
    if ~isfield(WaveformInfo, 'MaxChannel')
        WaveformInfo.MaxChannel = zeros(1, numel(UniqueIDConversion.UniqueID));
    end
end

function out = readPythonPickle(path)
%READPYTHONPICKLE Read a Python pickle via a temporary Python bridge.
    if ~exist(path, 'file')
        error('Pickle does not exist: %s', path);
    end
    tmpScript = fullfile(tempdir, 'read_python_pickle.py');
    fid = fopen(tmpScript, 'w');
    fprintf(fid, 'import pickle, json, sys\n');
    fprintf(fid, 'try:\n');
    fprintf(fid, '    import numpy as np\n');
    fprintf(fid, '    has_numpy = True\n');
    fprintf(fid, 'except ImportError:\n');
    fprintf(fid, '    has_numpy = False\n');
    fprintf(fid, '\n');
    fprintf(fid, 'class SafeEncoder(json.JSONEncoder):\n');
    fprintf(fid, '    def default(self, obj):\n');
    fprintf(fid, '        if has_numpy:\n');
    fprintf(fid, '            if isinstance(obj, np.ndarray):\n');
    fprintf(fid, '                return obj.tolist()\n');
    fprintf(fid, '            if isinstance(obj, (np.integer,)):\n');
    fprintf(fid, '                return int(obj)\n');
    fprintf(fid, '            if isinstance(obj, (np.floating,)):\n');
    fprintf(fid, '                return float(obj)\n');
    fprintf(fid, '            if isinstance(obj, (np.bool_,)):\n');
    fprintf(fid, '                return bool(obj)\n');
    fprintf(fid, '        if hasattr(obj, "__dict__"):\n');
    fprintf(fid, '            return obj.__dict__\n');
    fprintf(fid, '        if isinstance(obj, (set, frozenset)):\n');
    fprintf(fid, '            return list(obj)\n');
    fprintf(fid, '        try:\n');
    fprintf(fid, '            return str(obj)\n');
    fprintf(fid, '        except Exception:\n');
    fprintf(fid, '            return None\n');
    fprintf(fid, '\n');
    fprintf(fid, 'with open(sys.argv[1], "rb") as fh:\n');
    fprintf(fid, '    obj = pickle.load(fh)\n');
    fprintf(fid, 'if hasattr(obj, "__dict__"):\n');
    fprintf(fid, '    obj = obj.__dict__\n');
    fprintf(fid, 'sys.stdout.write(json.dumps(obj, cls=SafeEncoder))\n');
    fclose(fid);

    % Try python, then python3
    [status, outText] = system(sprintf('python "%s" "%s"', tmpScript, path));
    if status ~= 0
        [status, outText] = system(sprintf('python3 "%s" "%s"', tmpScript, path));
    end
    if status ~= 0
        error('Failed to read Python pickle: %s', strtrim(outText));
    end
    out = jsondecode(outText);
end

function out = readPythonNpz(path)
%READPYTHONNPZ Read a Python .npz file via a temporary Python bridge.
    if ~exist(path, 'file')
        error('NPZ does not exist: %s', path);
    end
    tmpScript = fullfile(tempdir, 'read_python_npz.py');
    fid = fopen(tmpScript, 'w');
    fprintf(fid, 'import numpy as np, json, sys\n');
    fprintf(fid, 'data = np.load(sys.argv[1], allow_pickle=True)\n');
    fprintf(fid, 'payload = {}\n');
    fprintf(fid, 'for key in data.files:\n');
    fprintf(fid, '    payload[key] = data[key].tolist()\n');
    fprintf(fid, 'print(json.dumps(payload))\n');
    fclose(fid);
    [status, outText] = system(sprintf('python "%s" "%s"', tmpScript, path));
    if status ~= 0
        error('Failed to read Python NPZ file.');
    end
    out = jsondecode(outText);
end
