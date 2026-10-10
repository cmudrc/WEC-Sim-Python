function most_short_source_baseline(end_time, output_step)
% Run the pinned MOST case with active input paths and selected output spacing.
if nargin == 0
    end_time = 10;
end
if nargin < 2
    output_step = 0.01;
end
assert(isfinite(end_time) && end_time > 0);
assert(isfinite(output_step) && output_step > 0);
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
contents = strrep(contents, old_end, sprintf('simu.endTime = %.15g;', end_time));
contents = strrep(contents, old_explorer, 'simu.explorer=''off'';');
contents = strrep(contents, old_output_step, ...
    sprintf('simu.dtOut = %.15g;', output_step));
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
body_acceleration = output.bodies(1).acceleration;
body_force_total = output.bodies(1).forceTotal;
body_force_excitation = output.bodies(1).forceExcitation;
body_force_radiation = output.bodies(1).forceRadiationDamping;
body_force_added_mass = output.bodies(1).forceAddedMass;
body_force_restoring = output.bodies(1).forceRestoring;
body_force_viscous = output.bodies(1).forceMorisonAndViscous;
body_force_linear_damping = output.bodies(1).forceLinearDamping;
mooring_position = output.mooring(1).position;
mooring_velocity = output.mooring(1).velocity;
mooring_force = output.mooring(1).forceMooring;
wave_time = output.wave.time;
wave_elevation = output.wave.elevation;
wave_omega = waves.omega;
wave_amplitude = waves.amplitude;
wave_phase = waves.phase;
wave_d_omega = waves.dOmega;
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
assert(all(isfinite(body_acceleration), 'all'));
assert(all(isfinite(body_force_total), 'all'));
assert(all(isfinite(mooring_force), 'all'));
assert(all(isfinite(rotor_speed), 'all'));
assert(all(isfinite(blade_pitch), 'all'));
assert(all(isfinite(generator_torque), 'all'));
assert(all(isfinite(azimuth), 'all'));
assert(all(isfinite(blade_aero_load), 'all'));
save(fullfile(root, 'matlab-most-short.mat'), ...
    'body_time', 'body_position', 'body_velocity', 'turbine_time', ...
    'rotor_speed', 'blade_pitch', 'generator_torque', 'wind_speed', ...
    'azimuth', 'blade_aero_load', 'tower_base_load', ...
    'body_acceleration', 'body_force_total', 'body_force_excitation', ...
    'body_force_radiation', 'body_force_added_mass', ...
    'body_force_restoring', 'body_force_viscous', ...
    'body_force_linear_damping', ...
    'mooring_position', 'mooring_velocity', 'mooring_force', ...
    'wave_time', 'wave_elevation', 'wave_omega', 'wave_amplitude', ...
    'wave_phase', 'wave_d_omega', '-v7');
end
