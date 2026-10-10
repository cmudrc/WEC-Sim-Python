function most_constant_wind_source_baseline
% Exercise the pinned MOST platform and active above-rated pitch control.
root = pwd;
input_file = fullfile(root, 'applications', 'MOST', 'wecSimInputFile.m');
contents = fileread(input_file);
old_flag = 'wind.constantWindFlag = 0;';
old_speed = 'windSpeed0=8;';
assert(contains(contents, old_flag) && contains(contents, old_speed));
contents = strrep(contents, old_flag, sprintf([ ...
    'wind.constantWindFlag = 1;\n' ...
    'wind.modules = [12 12];\n' ...
    'wind.timeBreakpoints = [0 10];\n' ...
    'wind.dt = 0.01;']));
contents = strrep(contents, old_speed, 'windSpeed0=12;');
fid = fopen(input_file, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', contents);
fclose(fid);

most_short_source_baseline(10);
movefile(fullfile(root, 'matlab-most-short.mat'), ...
    fullfile(root, 'matlab-most-constant-wind.mat'));
end
