function varargout=zmat(varargin)
%
% output=zmat(input)
%    or
% [output, info]=zmat(input, iscompress, method)
% output=zmat(input, info)
%
% A portable data compression/decompression toolbox for MATLAB/GNU Octave
% 
% author: Qianqian Fang <q.fang at neu.edu>
% initial version created on 04/30/2019
%
% input:
%      input: a char, non-complex numeric or logical vector or array
%      iscompress: (optional) if iscompress is 1, zmat compresses/encodes the input, 
%             if 0, it decompresses/decodes the input. Default value is 1.
%
%             if iscompress is set to a negative integer, (-iscompress) specifies
%             the compression level. For zlib/gzip, default level is 6 (1-9); for 
%             lzma/lzip, default level is 5 (1-9); for lz4hc, default level is 8 (1-16).
%             the default compression level is used if iscompress is set to 1.
%
%             zmat removes the trailing newline when iscompress=2 and method='base64'
%             all newlines are removed when iscompress=3 and method='base64'
%
%             if one defines iscompress as the info struct (2nd output of zmat), zmat 
%             will perform a decoding/decompression operation and recover the original
%             input using the info stored in the info structure.
%      method: (optional) compression method, currently, zmat supports the below methods
%             'zlib': zlib/zip based data compression (default)
%             'gzip': gzip formatted data compression
%             'lzip': lzip formatted data compression
%             'lzma': lzma formatted data compression
%             'lz4':  lz4 formatted data compression
%             'lz4hc':lz4hc (LZ4 with high-compression ratio) formatted data compression
%             'base64': encode or decode use base64 format
%
% output:
%      output: a uint8 row vector, storing the compressed or decompressed data; 
%             empty when an error is encountered
%      info: (optional) a struct storing additional info regarding the input data, may have
%            'type': the class of the input array
%            'size': the dimensions of the input array
%            'byte': the number of bytes per element in the input array
%            'method': a copy of the 3rd input indicating the encoding method
%            'status': the zlib/lzma/lz4 compression/decompression function return value, 
%                    including potential error codes; see documentation of the respective 
%                    libraries for details
%            'level': a copy of the iscompress flag; if non-zero, specifying compression 
%                    level, see above
%
% example:
%
%   [ss, info]=zmat(eye(5))
%   orig=zmat(ss,0)
%   orig=zmat(ss,info)
%   ss=char(zmat('zmat test',1,'base64'))
%   orig=char(zmat(ss,0,'base64'))
%
% -- this function is part of the ZMAT toolbox (http://github.com/fangq/zmat)
%

if(exist('zipmat')~=3 && exist('zipmat')~=2)
    error('zipmat mex file is not found. you must download the mex file or recompile');
end

if(nargin==0)
    fprintf(1,'Usage:\n\t[output,info]=zmat(input,iscompress,method);\nPlease run "help zmat" for more details.\n');
    return;
end

input=varargin{1};
iscompress=1;
zipmethod='zlib';

if(~(ischar(input) || islogical(input) || (isnumeric(input) && isreal(input))))
    error('input must be a char, non-complex numeric or logical vector or N-D array');
end

if(ischar(input))
    input=uint8(input);
end

if(nargin>1)
    iscompress=varargin{2};
    if(isstruct(varargin{2}))
        inputinfo=varargin{2};
        iscompress=0;
        zipmethod=inputinfo.method;
    end
end

if(nargin>2)
    zipmethod=varargin{3};
end

iscompress=round(iscompress);

if((strcmp(zipmethod,'zlib') || strcmp(zipmethod,'gzip')) && iscompress<=-10)
    iscompress=-9;
end

[varargout{1:max(1,nargout)}]=zipmat(input,iscompress,zipmethod);

if(strcmp(zipmethod,'base64') && iscompress>1)
    varargout{1}=char(varargout{1});
    if(iscompress==2)
        varargout{1}=regexprep(varargout{1},'\n$','');
    elseif(iscompress>2)
        varargout{1}=regexprep(varargout{1},'\n','');
    end
end

if(exist('inputinfo','var') && isfield(inputinfo,'type'))
        varargout{1}=typecast(varargout{1},inputinfo.type);
        varargout{1}=reshape(varargout{1},inputinfo.size);
end
