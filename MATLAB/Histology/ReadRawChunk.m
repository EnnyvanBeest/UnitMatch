function [data, fs, nChan] = ReadRawChunk(rawFile, sampleStart, nSamples)
% [data, fs, nChan] = ReadRawChunk(rawFile, sampleStart, nSamples)
% Reads nSamples samples (all channels) from a SpikeGLX .bin or mtscomp-compressed .cbin file,
% starting at sample sampleStart (0-based). For .cbin only the chunks that contain the range are
% read and decompressed (uses zmat, as bombcell's extractCbinData).
% data: nSamples x nChan int16; fs: sampling rate (Hz)

[fold,name,ext] = fileparts(rawFile);
if strcmpi(ext,'.cbin')
    ch = jsondecode(fileread(fullfile(fold,[name '.ch'])));
    nChan = ch.n_channels; fs = ch.sample_rate;
    bounds = ch.chunk_bounds(:); offsets = ch.chunk_offsets(:);
    sampleEnd = min(sampleStart + nSamples, bounds(end));
    c1 = find(bounds <= sampleStart,1,'last'); c2 = find(bounds >= sampleEnd,1,'first') - 1;
    info = struct('type',ch.dtype,'byte',2,'method',ch.algorithm,'status',1,'level',ch.comp_level);
    parts = cell(c2-c1+1,1);
    fid = fopen(rawFile,'r'); cleanup = onCleanup(@() fclose(fid));
    for c = c1:c2
        n = bounds(c+1) - bounds(c);
        fseek(fid,offsets(c),'bof');
        comp = fread(fid,offsets(c+1)-offsets(c),'*uint8');
        info.size = [n*nChan 1];
        d = reshape(zmat(comp,info),n,nChan);
        parts{c-c1+1} = cumsum(d,1); % time-differences per chunk, as in extractCbinData (stays int16)
    end
    data = cat(1,parts{:});
    first = sampleStart - bounds(c1) + 1;
    data = data(first:first+(sampleEnd-sampleStart)-1,:);
else
    % the meta file of this stream (ap or lf): the folder holds both, with different sampling rates
    [~,~,stream] = fileparts(name); stream = stream(2:end); % e.g. 'lf' for x.imec0.lf.bin
    if ismember(stream,{'ap','lf'}), meta = ReadMeta2(fold,stream); else, meta = ReadMeta2(fold); end
    nChan = str2double(meta.nSavedChans);
    if isfield(meta,'imSampRate'), fs = str2double(meta.imSampRate); else, fs = str2double(meta.niSampRate); end
    fid = fopen(rawFile,'r'); cleanup = onCleanup(@() fclose(fid));
    fseek(fid,double(sampleStart)*nChan*2,'bof');
    data = fread(fid,[nChan nSamples],'*int16')';
end
end
