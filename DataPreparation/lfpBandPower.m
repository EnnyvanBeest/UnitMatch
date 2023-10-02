

function [lfpByChannel, allPowerEst, F, allPowerVar] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand)
% function [lfpByChannel, allPowerEst, F] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand)
% Computes the power in particular bands, and across all frequencies, across the recording
% samples 10 segments of 10 sec each to compute these things. 

if ~isempty(freqBand) && ~iscell(freqBand)
    freqBand = {freqBand};
end
nF = length(freqBand);

nClips = 10;
clipDur = 10; % seconds

% load nClips one-sec samples
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;
sampStarts = round(linspace(lfpFs*10, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(lfpFs*clipDur);

mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});

allPowerEstByBand = zeros(nClips, nChansInFile, nF);
for n = 1:nClips
    fprintf(1, 'clip%d\n', n);
    % pull out the data
    thisDat = double(mmf.Data.x(:, (1:nClipSamps)+sampStarts(n)));
    
    % median subtract? 
%     thisDat = bsxfun(@minus, thisDat, median(thisDat));
    thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));
    
    [Pxx, F] = myTimePowerSpectrumMat(thisDat', lfpFs);
    
    if n==1
        allPowerEst = zeros(nClips, size(Pxx,1), size(Pxx,2));
    end
    allPowerEst(n,:,:) = Pxx;
        
    for f = 1:nF
        
        inclF = F>freqBand{f}(1) & F<=freqBand{f}(2);
        allPowerEstByBand(n,:, f) = mean(Pxx(inclF,:));
        
    end
end

if nF>0
    lfpByChannel = squeeze(mean(allPowerEstByBand, 1)); % mean across clips
else
    lfpByChannel = [];
end
allPowerVar = squeeze(var(allPowerEst,1));
allPowerEst = squeeze(mean(allPowerEst, 1));


function [Pxx, F] = myTimePowerSpectrumMat(x, Fs)
L = size(x,1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(x,[],[],NFFT,Fs);