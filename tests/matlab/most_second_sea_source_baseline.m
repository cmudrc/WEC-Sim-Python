function most_second_sea_source_baseline
% Exercise the published MOST configuration with a second irregular sea.
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

most_short_source_baseline;
movefile(fullfile(root, 'matlab-most-short.mat'), ...
    fullfile(root, 'matlab-most-sea6-seed2.mat'));
end
