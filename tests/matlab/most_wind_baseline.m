function most_wind_baseline
% Decode the checked-in MOST wind field with the published MATLAB reader.
source = fullfile(pwd, 'applications', 'MOST', 'mostData', ...
    'turbSim', 'WIND_8mps.bts');
addpath(fileparts(source));
[velocity, y, z, dt, hub_height] = readfile_BTS(source);
[nt, ~, ny, nz] = size(velocity);
selected_y = [1 6 12];
selected_z = [1 7 12];
traces = velocity(:, :, selected_y, selected_z);
spatial_sum = sum(sum(velocity, 4), 3);
spatial_square_sum = sum(sum(velocity.^2, 4), 3);
metadata = [nt ny nz dt hub_height];
save('matlab-most-wind.mat', 'metadata', 'y', 'z', 'traces', ...
    'spatial_sum', 'spatial_square_sum', '-v7');
end
