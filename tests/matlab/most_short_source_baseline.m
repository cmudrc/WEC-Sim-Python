function most_short_source_baseline
% Run the pinned MOST case through its first 10 s with active input paths.
root = pwd;
case_dir = fullfile(root, 'applications', 'MOST');
wind_grid = load(fullfile(root, 'matlab-most-advection.mat'));
SpatialDiscrUVW = wind_grid.values;
t = wind_grid.time;
Xdiscr = wind_grid.x;
Ydiscr = wind_grid.y;
Zdiscr = wind_grid.z;
save(fullfile(case_dir, 'mostData', 'turbSim', 'WIND_8mps.mat'), ...
    'SpatialDiscrUVW', 't', 'Xdiscr', 'Ydiscr', 'Zdiscr', '-v7.3');

input_file = fullfile(case_dir, 'wecSimInputFile.m');
contents = fileread(input_file);
old_end = 'simu.endTime = 1000;';
old_explorer = 'simu.explorer=''on'';';
old_output_step = 'simu.dtOut = 0.1;';
assert(contains(contents, old_end) && contains(contents, old_explorer) ...
    && contains(contents, old_output_step));
contents = strrep(contents, old_end, 'simu.endTime = 10;');
contents = strrep(contents, old_explorer, 'simu.explorer=''off'';');
contents = strrep(contents, old_output_step, 'simu.dtOut = 0.01;');
fid = fopen(input_file, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', contents);
fclose(fid);

cd(case_dir);
restore_dir = onCleanup(@() cd(root));
wecSim;
assert(numel(output.bodies) == 1 && numel(output.windTurbine) == 1);
body_time = output.bodies(1).time;
body_position = output.bodies(1).position;
body_velocity = output.bodies(1).velocity;
turbine_time = output.windTurbine(1).time;
rotor_speed = output.windTurbine(1).rotorSpeed;
blade_pitch = output.windTurbine(1).bladePitch;
generator_torque = output.windTurbine(1).genTorque;
wind_speed = output.windTurbine(1).windSpeed;
azimuth = output.windTurbine(1).azimuth;
blade_aero_load = cat(3, output.windTurbine(1).blade1AeroLoad, ...
    output.windTurbine(1).blade2AeroLoad, ...
    output.windTurbine(1).blade3AeroLoad);
tower_base_load = output.windTurbine(1).towerBaseLoad;
assert(all(isfinite(body_position), 'all'));
assert(all(isfinite(body_velocity), 'all'));
assert(all(isfinite(rotor_speed), 'all'));
assert(all(isfinite(blade_pitch), 'all'));
assert(all(isfinite(generator_torque), 'all'));
assert(all(isfinite(azimuth), 'all'));
assert(all(isfinite(blade_aero_load), 'all'));
save(fullfile(root, 'matlab-most-short.mat'), ...
    'body_time', 'body_position', 'body_velocity', 'turbine_time', ...
    'rotor_speed', 'blade_pitch', 'generator_torque', 'wind_speed', ...
    'azimuth', 'blade_aero_load', 'tower_base_load', '-v7');
end
