
mkdir(fullfile(savePath))

if ~exist(savePath,'dir')
    mkdir(fullfile(savePath))
end
disp([newline, 'saving quality metrics to ', savePath])
save(fullfile(savePath, 'qMetric.mat'), 'qMetric', '-v7.3')


qMetricSummary = table('Size',[length(qMetric.clusterID), 10],'VariableTypes',...
    {'double', 'double', 'double', 'double', 'double', 'double', 'double',...
    'double','double','double'},'VariableNames',...
    {'percSpikesMissing', 'clusterID', 'Fp', 'nSpikes', 'nPeaks', 'nTroughs', 'somatic', ...
    'waveformDuration', 'spatialDecaySlope', 'waveformBaseline'});
qMetricSummary.clusterID = qMetric.clusterID';
%qMetricSummary.percSpikesMissing = arrayfun(@(x) nanmean(qMetric.percSpikesMissing(qMetric.useTheseTimes{x})), 1:size(qMetric.percSpikesMissing,1));

qMetricSummary.percSpikesMissing = arrayfun(@(x) nanmean(qMetric.percSpikesMissing(x, qMetric.percSpikesMissing(x, :) <= param.maxPercSpikesMissing)), ...
    1:size(qMetric.percSpikesMissing, 1))';
qMetricSummary.Fp = arrayfun(@(x) nanmean(qMetric.Fp(x, qMetric.Fp(x, :) <= param.maxRPVviolations)), ...
    1:size(qMetric.percSpikesMissing, 1))';

qMetricSummary.percSpikesMissing = arrayfun( @(x) nanmean(qMetric.percSpikesMissing(x, qMetric.percSpikesMissing(x,:) <= param.maxPercSpikesMissing)), ...
    1:size(qMetric.percSpikesMissing,1))';
qMetricSummary.Fp = arrayfun( @(x) nanmean(qMetric.Fp(x, qMetric.Fp(x,:) <= param.maxRPVviolations)), ...
    1:size(qMetric.percSpikesMissing,1))';
qMetricSummary.nSpikes = qMetric.nSpikes';
qMetricSummary.nPeaks = qMetric.nPeaks';
qMetricSummary.nTroughs = qMetric.nTroughs';
qMetricSummary.somatic = qMetric.somatic';
qMetricSummary.waveformDuration = qMetric.waveformDuration';
qMetricSummary.spatialDecaySlope = qMetric.spatialDecaySlope';
qMetricSummary.waveformBaseline = qMetric.waveformBaseline';


parquetwrite([savePath, filesep, 'templates._jf_qMetrics.parquet'], qMetricSummary)

