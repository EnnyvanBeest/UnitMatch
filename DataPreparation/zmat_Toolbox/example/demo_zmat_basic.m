% addpath('../')

% compression
[dzip,info]=zmat(uint8(eye(5,5)))

% decompression
orig=reshape(zmat(dzip,0),info.size)

% base64 encoding and decoding
base64=zmat('zmat toolbox',1,'base64');
char(base64)

orig=zmat(base64,0,'base64');
char(orig)

% encode ND numeric array
orig=single(rand(5));
[Aencoded,info]=zmat(orig,1,'lzma')

% decode compressed ND array and restore size/type
Adecoded=zmat(Aencoded,info)

all(all(Adecoded==orig))
class(Adecoded)

