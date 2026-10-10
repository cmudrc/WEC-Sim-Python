function most_developed_sea_source_baseline(end_time)
% Run the pinned 6 m, seed-2 MOST sea beyond the wave ramp.
if nargin == 0
    end_time = 30;
end
root = pwd;
input_file = fullfile(root, 'applications', 'MOST', 'wecSimInputFile.m');
contents = fileread(input_file);
old_seed = 'waves.phaseSeed = 1;';
old_height = 'waves.height = 4;';
assert(contains(contents, old_seed) && contains(contents, old_height));
contents = strrep(contents, old_seed, 'waves.phaseSeed = 2;');
contents = strrep(contents, old_height, 'waves.height = 6;');
fid = fopen(input_file, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', contents);
fclose(fid);

most_short_source_baseline(end_time);
movefile(fullfile(root, 'matlab-most-short.mat'), ...
    fullfile(root, 'matlab-most-developed-sea.mat'));
end
