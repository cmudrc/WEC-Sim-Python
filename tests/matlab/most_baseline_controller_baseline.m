function most_baseline_controller_baseline
% Simulate the active published baseline controller block with generated data.
root = pwd;
data_dir = fullfile(root, 'applications', 'MOST', 'mostData', ...
    'windTurbine', 'control');
source = load(fullfile(data_dir, 'Control_IEA15MW.mat'), 'Ctrl');
steady = load(fullfile(data_dir, 'SteadyStates_IEA15MW.mat'), 'SteadyStates');
ss = steady.SteadyStates.ROSCO.SS;
omega0 = interp1(ss.WINDSPEED, ss.ROTSPD, 8);
pitch0 = interp1(ss.WINDSPEED, ss.BLADEPITCH, 8);
windTurbine.Baseline = source.Ctrl.Baseline;
windTurbine.omega0 = omega0;
windTurbine.bladepitch0 = pitch0;
assignin('base', 'windTurbine', windTurbine);

time = (0:0.01:40)';
speed = omega0 + 0.10*sin(0.7*time) ...
    + 0.22*(time >= 15) - 0.15*(time >= 30);
assignin('base', 'speed_signal', [time speed]);

library = fullfile(root, 'matlab-ref', 'source', 'lib', 'MOST', 'MOST_Lib.slx');
load_system(library);
model = 'most_baseline_controller_harness';
new_system(model);
cleanup = onCleanup(@() close_system(model, 0));
add_block('simulink/Sources/From Workspace', [model '/RotorSpeed'], ...
    'VariableName', 'speed_signal');
add_block(['MOST_Lib' ...
    '/Wind turbine/Aerodynamics + Control/Control/Baseline'], ...
    [model '/Baseline']);
set_param([model '/Baseline'], 'VariantControl', '');
add_block('simulink/Sinks/To Workspace', [model '/GenTorque'], ...
    'VariableName', 'genTorque', 'SaveFormat', 'Array');
add_block('simulink/Sinks/To Workspace', [model '/BladePitch'], ...
    'VariableName', 'bladePitch', 'SaveFormat', 'Array');
add_line(model, 'RotorSpeed/1', 'Baseline/1');
add_line(model, 'Baseline/1', 'GenTorque/1');
add_line(model, 'Baseline/2', 'BladePitch/1');
set_param(model, 'Solver', 'ode4', 'FixedStep', '0.01', ...
    'StartTime', '0', 'StopTime', '40', 'SaveTime', 'on');
output = sim(model, 'ReturnWorkspaceOutputs', 'on');
torque = output.get('genTorque');
pitch = output.get('bladePitch');
assert(numel(torque) == numel(time));
assert(numel(pitch) == numel(time));
assert(all(isfinite(torque)) && all(isfinite(pitch)));
control = source.Ctrl.Baseline;
save(fullfile(root, 'matlab-most-controller.mat'), ...
    'control', 'omega0', 'pitch0', 'time', 'speed', 'torque', 'pitch', '-v7');
end
