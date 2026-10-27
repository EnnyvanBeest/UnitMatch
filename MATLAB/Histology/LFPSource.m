function [rawFile, Shank, Depth, apToRaw] = LFPSource(apFile)
% [rawFile, Shank, Depth, apToRaw] = LFPSource(apFile)
% The file to read LFP from for a SpikeGLX AP file, and the position of its channels.
% NP2.0 records wideband: the AP file itself. NP1.0 (AP high-passed): the .lf.bin/.lf.cbin next to
% it (rawFile empty if there is none).
% Shank (0-based; shanks 250 um apart) and Depth (um) per channel, in raw channel order without
% the sync channel. apToRaw: factor from AP samples to rawFile samples (1, or 1/12 for NP1 LF).
[fold,name,~] = fileparts(apFile);
meta = ReadMeta2(fold,'ap');
rawFile = apFile; apToRaw = 1;
if isfield(meta,'imDatPrb_type') && ismember(str2double(meta.imDatPrb_type),[0 1020 1030 1100 1120 1121 1122 1123 1200 1300]) % NP1-like: separate LF file
    lf = [dir(fullfile(fold,[strrep(name,'.ap','.lf') '.cbin'])); dir(fullfile(fold,[strrep(name,'.ap','.lf') '.bin']))];
    if isempty(lf)
        warning('LFPSource:NoLF','No LF file for %s',apFile); rawFile = '';
    else
        rawFile = fullfile(lf(1).folder,lf(1).name);
        lfMeta = ReadMeta2(lf(1).folder,'lf');
        apToRaw = str2double(lfMeta.imSampRate)/str2double(meta.imSampRate);
    end
end
chanPos = ChannelIMROConversion(fold,0,0);
Shank = floor(chanPos(:,1)/250); Depth = chanPos(:,2);
end
