function tint = loadintants(fname,fs)

if nargin < 2
    fs = 30000;
end

dig = readNPY(fname);
dig = int8(dig(1,:));
digd = diff(dig);
ixa = digd==1;
ixb = digd==-1;
dt = 1/fs;
tint = 0:dt:dt*length(dig);
tint = tint(ixa|ixb);