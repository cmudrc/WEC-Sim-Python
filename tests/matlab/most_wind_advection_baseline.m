function most_wind_advection_baseline(frames)
% Run the published MOST wind preprocessor on a bounded BTS excerpt.
if nargin == 0
    frames = 1000;
end
root = pwd;
source_dir = fullfile(root, 'applications', 'MOST', 'mostData', 'turbSim');
output_dir = fullfile(root, 'matlab-most-advection');
mkdir(output_dir);
copyfile(fullfile(source_dir, 'RunTurbsim.m'), output_dir);
copyfile(fullfile(source_dir, 'readfile_BTS.m'), output_dir);

input_file = fullfile(source_dir, 'WIND_8mps.bts');
fid = fopen(input_file, 'rb');
assert(fid > 0);
header = fread(fid, 70, '*uint8');
assert(numel(header) == 70);
nchar = double(typecast(header(67:70), 'int32'));
nz = double(typecast(header(3:6), 'int32'));
ny = double(typecast(header(7:10), 'int32'));
ntwr = double(typecast(header(11:14), 'int32'));
body = fread(fid, nchar + frames*2*3*(ny*nz + ntwr), '*uint8');
fclose(fid);
assert(numel(body) == nchar + frames*2*3*(ny*nz + ntwr));
header(15:18) = reshape(typecast(int32(frames), 'uint8'), [], 1);
fid = fopen(fullfile(output_dir, 'WIND_8mps.bts'), 'wb');
fwrite(fid, header, 'uint8');
fwrite(fid, body, 'uint8');
fclose(fid);

addpath(output_dir);
cd(output_dir);
restore_dir = onCleanup(@() cd(root));
Wind = RunTurbsim();
values = Wind.SpatialDiscrUVW;
time = Wind.t;
x = Wind.Xdiscr;
y = Wind.Ydiscr;
z = Wind.Zdiscr;
save(fullfile(root, 'matlab-most-advection.mat'), ...
    'values', 'time', 'x', 'y', 'z', '-v7');
end
