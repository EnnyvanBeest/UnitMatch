% =============================================================
% Simple helper functions and MATLAB structures demonstrating
% how to read and manipulate SpikeGLX meta and binary files.
%
% The most important part of the demo is ReadMeta().
% Please read the comments for that function. Use of
% the 'meta' structure will make your data handling
% much easier!
%
function [lfpByChannel, allPowerEst, F, allPowerVar] = lfpBandPowerNP2(lfpFilename,freqBand)
% function [lfpByChannel, allPowerEst, F] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand)
% Computes the power in particular bands (freqBand), and across all frequencies, across the recording
% samples 10 segments of 10 sec each to compute these things. 
%Downsample function
DownsampleAvg = @(y,n) nanmean(reshape([y(:); nan(mod(-numel(y),n),1)],n,[]));
compressed = 0;
if any(strfind(lfpFilename,'.cbin'))
    compressed = 1;
end
if nargin==0 || ~exist(lfpFilename)
    [binName,path] = uigetfile('*.bin', 'Select Binary File');
else
    
    [path,binName] = fileparts(lfpFilename);
    if ~compressed && isempty(strfind(binName,'.bin')) %somehow it cuts the bin off
        binName =[binName '.bin'];
    elseif compressed &&  isempty(strfind(binName,'.cbin')) 
        binName =[binName '.cbin'];
    end
end

if ~isempty(freqBand) && ~iscell(freqBand)
    freqBand = {freqBand};
end
nF = length(freqBand);

nClips = 5;
clipDur = 10; % seconds

% Parse the corresponding metafile
meta = ReadMeta(binName, path);
nSamps = str2num(meta.fileTimeSecs)*str2num(meta.imSampRate);
sampStarts = round(linspace(str2num(meta.imSampRate)*10, nSamps, nClips+1)); % skip first 10 secs
nClipSamps = round(str2num(meta.imSampRate)*clipDur);
nChansInFile = str2num(meta.nSavedChans);
dataType = 'A';         %set to 'A' for analog, 'D' for digital data

% For an analog channel: gain correct saved channel ch (1-based for MATLAB).
ch = 1;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;

% Read these lines in dw (0-based).
% For 3B2 imec data: the sync pulse is stored in line 6.
% May be 1 or more line indices.
dLineList = [0,1,6];

% Get first 30 seconds of data
allPowerEstByBand = zeros(nClips, nChansInFile, nF);
FreqNeeded = 1000;
for n = 1:nClips
    fprintf(1, 'clip%d\n', n);
    if  ~compressed
        dataArray = ReadBin(sampStarts(n), nClipSamps, meta, binName, path); %This is stored in int16
    else
        dataArray = pyrunfile("Ephys_Reader_FromMatlab.py","chunk",...
            datapath = strrep(fullfile(path,binName),'\','/'),start_time=sampStarts(n),end_time=sampStarts(n)+nClipSamps); %0-indexed!!
        dataArray=uint16(dataArray);


    end
    % convert to double
    dataArray = double(dataArray);
    if dataType == 'A'
        if strcmp(meta.typeThis, 'imec')
            dataArray = GainCorrectIM(dataArray, [ch], meta);
        else
            dataArray = GainCorrectNI(dataArray, [ch], meta);
        end
%         plot(dataArray(ch,:));
    else
        digArray = ExtractDigital(dataArray, meta, dw, dLineList);
        for i = 1:numel(dLineList)
%             plot(digArray(i,:));
%             hold on
        end
%         hold off
    end
    dataArray = arrayfun(@(ch) DownsampleAvg(dataArray(ch,:),SampRate(meta)./FreqNeeded),1:nChansInFile,'UniformOutput',0);
    dataArray = cat(1,dataArray{:});
    
    %Mean subtract
    dataArray = bsxfun(@minus,dataArray,mean(dataArray,2));
    [Pxx, F] = myTimePowerSpectrumMat(dataArray', FreqNeeded);
    
    if n==1
        allPowerEst = zeros(nClips, size(Pxx,1), size(Pxx,2));
    end
    allPowerEst(n,:,:) = Pxx;
    for f = 1:nF
        
        inclF = F>freqBand{f}(1) & F<=freqBand{f}(2);
        allPowerEstByBand(n,:, f) = mean(Pxx(inclF,:));
        
    end
    clear dataArray
    clear Pxx
    
end
if nF>0
    lfpByChannel = squeeze(mean(allPowerEstByBand, 1)); % mean across clips
else
    lfpByChannel = [];
end
allPowerVar = squeeze(var(allPowerEst,1));
allPowerEst = squeeze(mean(allPowerEst, 1));

end % DemoReadSGLXData


% =========================
% General Utility Functions
% =========================


% =========================================================
% Parse ini file returning a structure whose field names
% are the metadata left-hand-side tags, and whose right-
% hand-side values are MATLAB strings. We remove any
% leading '~' characters from tags because MATLAB uses
% '~' as an operator.
%
% If you're unfamiliar with structures, the benefit
% is that after calling this function you can refer
% to metafile items by name. For example:
%
%   meta.fileCreateTime  // file create date and time
%   meta.nSavedChans     // channels per timepoint
%
% All of the values are MATLAB strings, but you can
% obtain a numeric value using str2double(meta.nSavedChans).
% More complicated parsing of values is demonstrated in the
% utility functions below.
%
function [meta] = ReadMeta(binName, path)

    % Create the matching metafile name
    [dumPath,name,dumExt] = fileparts(binName);
    metaName = strcat(name, '.meta');

    % Parse ini file into cell entries C{1}{i} = C{2}{i}
    fid = fopen(fullfile(path, metaName), 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
    C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
    fclose(fid);

    % New empty struct
    meta = struct();

    % Convert each cell entry into a struct entry
    for i = 1:length(C{1})
        tag = C{1}{i};
        if tag(1) == '~'
            % remake tag excluding first character
            tag = sprintf('%s', tag(2:end));
        end
        meta = setfield(meta, tag, C{2}{i});
    end
end % ReadMeta


% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [nChan,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.
%
function dataArray = ReadBin(samp0, nSamp, meta, binName, path)

    nChan = str2double(meta.nSavedChans);

    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    samp0 = max(samp0, 0);
    nSamp = min(nSamp, nFileSamp - samp0);

    sizeA = [nChan, nSamp];

    fid = fopen(fullfile(path, binName), 'rb');
    fseek(fid, samp0 * 2 * nChan, 'bof');
    dataArray = fread(fid, sizeA, 'int16=>uint16');
    fclose(fid);
end % ReadBin


% =========================================================
% Return an array [lines X timepoints] of uint8 values for
% a specified set of digital lines.
%
% - dwReq is the one-based index into the saved file of the
%    16-bit word that contains the digital lines of interest.
% - dLineList is a zero-based list of one or more lines/bits
%    to scan from word dwReq.
%
function digArray = ExtractDigital(dataArray, meta, dwReq, dLineList)
    % Get channel index of requested digital word dwReq
    if strcmp(meta.typeThis, 'imec')
        [AP, LF, SY] = ChannelCountsIM(meta);
        if SY == 0
            fprintf('No imec sync channel saved\n');
            digArray = [];
            return;
        else
            digCh = AP + LF + dwReq;
        end
    else
        [MN,MA,XA,DW] = ChannelCountsNI(meta);
        if dwReq > DW
            fprintf('Maximum digital word in file = %d\n', DW);
            digArray = [];
            return;
        else
            digCh = MN + MA + XA + dwReq;
        end
    end
    [~,nSamp] = size(dataArray);
    digArray = zeros(numel(dLineList), nSamp, 'uint8');
    for i = 1:numel(dLineList)
        digArray(i,:) = bitget(dataArray(digCh,:), dLineList(i)+1, 'int16');
    end
end % ExtractDigital


% =========================================================
% Return sample rate as double.
%
function srate = SampRate(meta)
    if strcmp(meta.typeThis, 'imec')
        srate = str2double(meta.imSampRate);
    else
        srate = str2double(meta.niSampRate);
    end
end % SampRate


% =========================================================
% Return a multiplicative factor for converting 16-bit
% file data to voltage. This does not take gain into
% account. The full conversion with gain is:
%
%   dataVolts = dataInt * fI2V / gain.
%
% Note that each channel may have its own gain.
%
function fI2V = Int2Volts(meta)
    if strcmp(meta.typeThis, 'imec')
        if isfield(meta,'imMaxInt')
            maxInt = str2num(meta.imMaxInt);
        else
            maxInt = 512;
        end
        fI2V = str2double(meta.imAiRangeMax) / maxInt;
    else
        fI2V = str2double(meta.niAiRangeMax) / 32768;
    end
end % Int2Volts


% =========================================================
% Return array of original channel IDs. As an example,
% suppose we want the imec gain for the ith channel stored
% in the binary data. A gain array can be obtained using
% ChanGainsIM() but we need an original channel index to
% do the look-up. Because you can selectively save channels
% the ith channel in the file isn't necessarily the ith
% acquired channel, so use this function to convert from
% ith stored to original index.
%
% Note: In SpikeGLX channels are 0-based, but MATLAB uses
% 1-based indexing, so we add 1 to the original IDs here.
%
function chans = OriginalChans(meta)
    if strcmp(meta.snsSaveChanSubset, 'all')
        chans = (1:str2double(meta.nSavedChans));
    else
        chans = str2num(meta.snsSaveChanSubset);
        chans = chans + 1;
    end
end % OriginalChans


% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
function [AP,LF,SY] = ChannelCountsIM(meta)
    M = str2num(meta.snsApLfSy);
    AP = M(1);
    LF = M(2);
    SY = M(3);
end % ChannelCountsIM

% =========================================================
% Return counts of each nidq channel type that compose
% the timepoints stored in binary file.
%
function [MN,MA,XA,DW] = ChannelCountsNI(meta)
    M = str2num(meta.snsMnMaXaDw);
    MN = M(1);
    MA = M(2);
    XA = M(3);
    DW = M(4);
end % ChannelCountsNI


% =========================================================
% Return gain for ith channel stored in the nidq file.
%
% ichan is a saved channel index, rather than an original
% (acquired) index.
%
function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
    if ichan <= savedMN
        gain = str2double(meta.niMNGain);
    elseif ichan <= savedMN + savedMA
        gain = str2double(meta.niMAGain);
    else
        gain = 1;
    end
end % ChanGainNI


% =========================================================
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
function [APgain,LFgain] = ChanGainsIM(meta)

    if isfield(meta,'imDatPrb_type')
        probeType = str2num(meta.imDatPrb_type);
    else
        probeType = 0;
    end
    if (probeType == 21) || (probeType == 24)
        [AP,LF,~] = ChannelCountsIM(meta);
        % NP 2.0; APgain = 80 for all channels
        APgain = zeros(AP,1,'double');
        APgain = APgain + 80;
        % No LF channels, set gain = 0
        LFgain = zeros(LF,1,'double');
    else
        % 3A or 3B data?
        % 3A metadata has field "typeEnabled" which was replaced
        % with "typeImEnabled" and "typeNiEnabled" in 3B.
        % The 3B imro table has an additional field for the
        % high pass filter enabled/disabled
        if isfield(meta,'typeEnabled')
            % 3A data
            C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
        else
            % 3B data
            C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
        end
        APgain = double(cell2mat(C(1)));
        LFgain = double(cell2mat(C(2)));
    end
end % ChanGainsIM


% =========================================================
% Having acquired a block of raw nidq data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:MN]      % all MN chans (MN from ChannelCountsNI)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectNI(dataArray, chanList, meta)

    [MN,MA] = ChannelCountsNI(meta);
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        conv = fI2V / ChanGainNI(j, MN, MA, meta);
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end


% =========================================================
% Having acquired a block of raw imec data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:AP]      % all AP chans (AP from ChannelCountsIM)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectIM(dataArray, chanList, meta)

    % Look up gain with acquired channel ID
    chans = OriginalChans(meta);
    [APgain,LFgain] = ChanGainsIM(meta);
    nAP = length(APgain);
    nNu = nAP * 2;

    % Common conversion factor
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        k = chans(j);       % acquisition index
        if k <= nAP
            conv = fI2V / APgain(k);
        elseif k <= nNu
            conv = fI2V / LFgain(k - nAP);
        else
            continue;
        end
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end

function [Pxx, F] = myTimePowerSpectrumMat(x, Fs)
L = size(x,1);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(x,[],[],NFFT,Fs);
end
