function WriteRaw(fname,Raw,type,varargin)

% Usage WriteRaw('image.raw',Raw,'int16');
% WriteRaw('image.raw',Raw, 'int16','BigEndian');


fid = fopen(fname, 'w');
i=1;

if strcmp(varargin,'BigEndian')
    mformat = 'ieee-be';
elseif strcmp(varargin,'LittleEndian')    
    mformat = 'ieee-le';
else 
    mformat = [];
end

zdim = size(Raw,3);
if isempty(mformat)
    while i<=zdim
        fwrite(fid, Raw(:,:,i), type);
        i=i+1;
    end
else
    while i<=zdim
        fwrite(fid, Raw(:,:,i), type, mformat);
        i=i+1;
    end
end

fclose(fid);