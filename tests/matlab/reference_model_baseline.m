function reference_model_baseline(model)
% Run the published RM3, OSWEC, or Sphere reference cases and save numeric output.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
outDir = fullfile(repoRoot, 'matlab-reference-model-output');
if ~isfolder(outDir)
    mkdir(outDir);
end
if any(string(model) == ["RM3_MCR_ARRAY", "RM3_MCR_EXCEL", "RM3_MCR_MAT"])
    run_rm3_mcr_baseline(repoRoot, outDir, string(model));
    return;
end
if string(model) == "RM3_MCR_SEASTATE"
    run_rm3_mcr_seastate_baseline(repoRoot, outDir);
    return;
end

switch string(model)
    case "RM3"
        cases = string(model);
        caseDirs = string(fullfile(repoRoot, 'matlab-ref', 'examples', model));
    case "OSWEC"
        cases = string(model);
        caseDirs = string(fullfile(repoRoot, 'matlab-ref', 'examples', model));
        inputFile = fullfile(caseDirs, 'wecSimInputFile.m');
        contents = fileread(inputFile);
        oldSetting = 'waves.spread = [0.1,0.2,0.7];';
        assert(contains(contents, oldSetting), 'The pinned OSWEC input changed');
        contents = strrep(contents, oldSetting, ...
            sprintf('%s\nwaves.phaseSeed = 1;', oldSetting));
        fid = fopen(inputFile, 'w');
        assert(fid ~= -1, 'Could not seed the pinned OSWEC input');
        fprintf(fid, '%s', contents);
        fclose(fid);
    case "Sphere"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', 'Sphere', 'hydroData');
        cd(hydroDir);
        if ~isfile('sphere.h5')
            bemio;
        end
        cases = ["0m", "1m", "1m-ME", "3m", "5m"];
        caseDirs = fullfile(repoRoot, 'applications', 'Free_Decay', cases);
    case "RM3_B2B"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', 'RM3', 'hydroData');
        cd(hydroDir);
        if ~isfile('rm3.h5')
            bemio;
        end
        cases = ["B2B_Case1", "B2B_Case2", "B2B_Case3", "B2B_Case4"];
        caseDirs = fullfile(repoRoot, 'applications', ...
            'Body-to-Body_Interactions', cases);
    case "RM3_END_STOPS"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'RM3', 'hydroData');
        cd(hydroDir);
        if ~isfile('rm3.h5')
            bemio;
        end
        cases = "End_Stops";
        caseDirs = string(fullfile(repoRoot, 'applications', 'End_Stops'));
    case {"RM3_END_STOPS_STEP", "RM3_END_STOPS_STEP_FINE"}
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'RM3', 'hydroData');
        cd(hydroDir);
        if ~isfile('rm3.h5')
            bemio;
        end
        sourceDir = fullfile(repoRoot, 'applications', 'End_Stops');
        if string(model) == "RM3_END_STOPS_STEP"
            cases = "End_Stops_dt005";
            stepSize = 0.05;
        else
            cases = "End_Stops_dt0025";
            stepSize = 0.025;
        end
        caseDirs = string(fullfile(repoRoot, 'applications', cases));
        [copied, copyMessage] = copyfile(sourceDir, caseDirs);
        assert(copied, copyMessage);
        inputFile = fullfile(caseDirs, 'wecSimInputFile.m');
        contents = fileread(inputFile);
        assert(contains(contents, 'simu.dt = 0.1;') && ...
            contains(contents, 'simu.endTime = 400;'), ...
            'The pinned End_Stops time settings changed');
        contents = strrep(contents, 'simu.dt = 0.1;', ...
            sprintf('simu.dt = %.3f;', stepSize));
        contents = strrep(contents, 'simu.endTime = 400;', ...
            'simu.endTime = 120;');
        fid = fopen(inputFile, 'w');
        assert(fid ~= -1, 'Could not write the End_Stops step input');
        fprintf(fid, '%s', contents);
        fclose(fid);
    case "OSWEC_Nonhydro"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'OSWEC', 'hydroData');
        cd(hydroDir);
        if ~isfile('oswec.h5')
            bemio;
        end
        cases = "Nonhydro";
        caseDirs = string(fullfile(repoRoot, 'applications', 'Nonhydro_Body'));
    case "RM3_PTO_Extension"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'RM3', 'hydroData');
        cd(hydroDir);
        if ~isfile('rm3.h5')
            bemio;
        end
        cases = ["float", "spar"];
        caseDirs = fullfile(repoRoot, 'applications', ...
            'RM3_PTO_Extension', cases);
    case "RM3_Radiation_Options"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'RM3', 'hydroData');
        cd(hydroDir);
        if ~isfile('rm3.h5')
            bemio;
        end
        cases = ["constant", "FIR", "state_space", "convolution"];
        caseDirs = repmat(string(fullfile(repoRoot, 'applications', ...
            'Radiation_Force_Options')), size(cases));
    case "RM3_MCR"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', 'RM3', 'hydroData');
        cd(hydroDir);
        if ~isfile('rm3.h5')
            bemio;
        end
        sourceDir = fullfile(repoRoot, 'applications', 'Multiple_Condition_Runs', 'RM3_MCROPT1');
        cases = strings(1, 8);
        caseDirs = strings(1, 8);
        iCase = 0;
        for damping = [1200000, 2400000]
            for period = [6, 8]
                for height = [1.5, 2.5]
                    iCase = iCase + 1;
                    cases(iCase) = sprintf('H%d_T%d_D%d', ...
                        round(10 * height), period, round(damping / 100000));
                    caseDirs(iCase) = fullfile(repoRoot, 'applications', ...
                        'Multiple_Condition_Runs', 'paired_' + cases(iCase));
                    [copied, copyMessage] = copyfile(sourceDir, caseDirs(iCase));
                    assert(copied, copyMessage);
                    inputFile = fullfile(caseDirs(iCase), 'wecSimInputFile.m');
                    contents = fileread(inputFile);
                    oldSettings = ["waves.height = 1.5:1:2.5;", ...
                        "waves.period = 6:2:8;", ...
                        "pto(1).damping=1200000:1200000:2400000;"];
                    newSettings = ["waves.height = " + string(height) + ";", ...
                        "waves.period = " + string(period) + ";", ...
                        "pto(1).damping=" + string(damping) + ";"];
                    for iSetting = 1:numel(oldSettings)
                        assert(contains(contents, oldSettings(iSetting)), ...
                            'The pinned RM3 MCR input changed');
                        contents = strrep(contents, char(oldSettings(iSetting)), ...
                            char(newSettings(iSetting)));
                    end
                    fid = fopen(inputFile, 'w');
                    assert(fid ~= -1, 'Could not write the RM3 MCR input');
                    fprintf(fid, '%s', contents);
                    fclose(fid);
                end
            end
        end
    case "Sphere_Passive"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', 'Sphere', 'hydroData');
        cd(hydroDir);
        if ~isfile('sphere.h5')
            bemio;
        end
        cases = "Passive_P";
        caseDirs = string(fullfile(repoRoot, 'applications', 'Controls', 'Passive (P)'));
    case "Sphere_PTO_Config"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', 'Sphere', 'hydroData');
        cd(hydroDir);
        if ~isfile('sphere.h5')
            bemio;
        end
        sourceDir = fullfile(repoRoot, 'applications', 'Controls', 'Passive (P)');
        caseDir = fullfile(repoRoot, 'applications', 'Controls', 'Passive_Config');
        [copied, copyMessage] = copyfile(sourceDir, caseDir);
        assert(copied, copyMessage);
        inputFile = fullfile(caseDir, 'wecSimInputFile.m');
        contents = fileread(inputFile);
        oldSettings = ["pto(1).stiffness = 0;", "pto(1).damping = 0;", ...
            "pto(1).location = [0 0 0];"];
        newSettings = ["pto(1).stiffness = 50000;", ...
            "pto(1).damping = 100000;", "pto(1).location = [1 0 0];"];
        for iSetting = 1:numel(oldSettings)
            assert(contains(contents, oldSettings(iSetting)), ...
                'The pinned passive-controller input changed');
            contents = strrep(contents, char(oldSettings(iSetting)), ...
                char(newSettings(iSetting)));
        end
        fid = fopen(inputFile, 'w');
        assert(fid ~= -1, 'Could not write the configured PTO input');
        fprintf(fid, '%s', contents);
        fclose(fid);
        cases = "Configured";
        caseDirs = string(caseDir);
    case "Sphere_Reactive_PI"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'Sphere', 'hydroData');
        cd(hydroDir);
        if ~isfile('sphere.h5')
            bemio;
        end
        cases = "Reactive_PI";
        caseDirs = string(fullfile(repoRoot, 'applications', 'Controls', ...
            'Reactive (PI)'));
    case "Sphere_Declutching"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'Sphere', 'hydroData');
        cd(hydroDir);
        if ~isfile('sphere.h5')
            bemio;
        end
        cases = "Declutching";
        caseDirs = string(fullfile(repoRoot, 'applications', 'Controls', ...
            'Declutching'));
    case "Sphere_Latching"
        hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
            'Sphere', 'hydroData');
        cd(hydroDir);
        if ~isfile('sphere.h5')
            bemio;
        end
        cases = "Latching";
        caseDirs = string(fullfile(repoRoot, 'applications', 'Controls', ...
            'Latching'));
    otherwise
        error('Unknown reference model: %s', model);
end

for iCase = 1:numel(cases)
    cd(caseDirs(iCase));
    if string(model) == "RM3_Radiation_Options"
        [output, simu] = run_radiation_option(cases(iCase));
    else
        wecSim;
    end
    assert(exist('output', 'var') == 1 && ~isempty(output.bodies), ...
        'The MATLAB case produced no body output');
    if any(string(model) == ["RM3_END_STOPS", "RM3_END_STOPS_STEP", ...
                             "RM3_END_STOPS_STEP_FINE"])
        stops = pto(1).hardStops;
        assert(strcmp(stops.upperLimitSpecify, 'on') && ...
            strcmp(stops.lowerLimitSpecify, 'on') && ...
            stops.upperLimitBound == 0.6 && stops.lowerLimitBound == -0.6 && ...
            stops.upperLimitStiffness == 1e8 && ...
            stops.lowerLimitStiffness == 1e8 && ...
            stops.upperLimitDamping == 0 && stops.lowerLimitDamping == 0 && ...
            stops.upperLimitTransitionRegionWidth == 1e-4 && ...
            stops.lowerLimitTransitionRegionWidth == 1e-4, ...
            'The pinned End_Stops effective hard-stop settings changed');
        writematrix([stops.lowerLimitBound, stops.upperLimitBound, ...
            stops.lowerLimitStiffness, stops.upperLimitStiffness, ...
            stops.lowerLimitDamping, stops.upperLimitDamping, ...
            stops.lowerLimitTransitionRegionWidth, ...
            stops.upperLimitTransitionRegionWidth], ...
            fullfile(outDir, string(model) + "_settings.csv"));
        if any(string(model) == ["RM3_END_STOPS_STEP", ...
                                 "RM3_END_STOPS_STEP_FINE"])
            assert(simu.dt == stepSize && simu.endTime == 120, ...
                'The End_Stops step diagnostic did not use its requested settings');
            writematrix([simu.dt, simu.endTime, simu.rampTime], ...
                fullfile(outDir, string(model) + "_run_settings.csv"));
        end
    end
    if string(model) == "OSWEC"
        % The published case shuffles its random phase. The baseline sets
        % phaseSeed=1 and saves the realized components for paired checks.
        assert(waves.phaseSeed == 1, 'OSWEC phase seed was not applied');
        components = [waves.omega(:), waves.amplitude(:), ...
            waves.dOmega(:), waves.phase];
        assert(all(isfinite(components), 'all'), 'Nonfinite wave components');
        writematrix(components, fullfile(outDir, 'OSWEC_wave_components.csv'));
        writematrix([waves.direction(:), waves.spread(:)], ...
            fullfile(outDir, 'OSWEC_wave_directions.csv'));
        writematrix(waves.waveAmpTime, fullfile(outDir, 'OSWEC_wave_elevation.csv'));
    end
    for iBody = 1:numel(output.bodies)
        response = output.bodies(iBody);
        values = [response.time(:), response.position, response.velocity, ...
            response.forceTotal, response.forceExcitation];
        if any(string(model) == ["RM3_END_STOPS", "RM3_END_STOPS_STEP", ...
                                 "RM3_END_STOPS_STEP_FINE"])
            forceValues = [response.time(:), response.forceRadiationDamping, ...
                response.forceAddedMass, response.forceRestoring, ...
                response.acceleration];
            assert(all(isfinite(forceValues), 'all'), ...
                'End-stop hydrodynamic forces contain nonfinite values');
            writematrix(forceValues, fullfile(outDir, sprintf( ...
                '%s_body%d_forces.csv', model, iBody)));
            hydroForce = body(iBody).hydroForce.hf1;
            assert(isfield(hydroForce.storage, 'hydroForce_fAddedMass'), ...
                'End-stop case is missing its applied added-mass matrix');
            writematrix(hydroForce.fAddedMass, fullfile(outDir, sprintf( ...
                '%s_body%d_added_mass_original.csv', model, iBody)));
            writematrix(hydroForce.storage.hydroForce_fAddedMass, ...
                fullfile(outDir, sprintf( ...
                '%s_body%d_added_mass_applied.csv', model, iBody)));
            properties = [body(iBody).mass, hydroForce.mass, ...
                hydroForce.adjustedMass, body(iBody).inertia, ...
                hydroForce.adjustedInertia];
            writematrix(properties, fullfile(outDir, sprintf( ...
                '%s_body%d_mass_properties.csv', model, iBody)));
        end
        if string(model) == "RM3_Radiation_Options"
            values = [values, response.forceRadiationDamping, ...
                response.forceAddedMass];
        end
        assert(all(isfinite(values), 'all'), 'The MATLAB response contains nonfinite values');
        filename = sprintf('%s_%s_body%d.csv', model, cases(iCase), iBody);
        writematrix(values, fullfile(outDir, filename));
    end
    if ismember(string(model), ["Sphere_Passive", "Sphere_PTO_Config", ...
            "Sphere_Reactive_PI", "Sphere_Declutching", "Sphere_Latching"])
        assert(exist('controller1_out', 'var') == 1, ...
            'The Sphere controller produced no logged output');
        controllerValues = [controller1_out.time(:), controller1_out.signals.values];
        assert(size(controllerValues, 2) == 13, ...
            'Expected six force and six power components from Sphere controller');
        assert(all(isfinite(controllerValues), 'all'), ...
            'The passive controller output contains nonfinite values');
        writematrix(controllerValues, fullfile(outDir, string(model) + "_controller.csv"));
    end
    if isstruct(output.ptos) && isfield(output.ptos, 'time')
        for iPto = 1:numel(output.ptos)
            response = output.ptos(iPto);
            values = [response.time(:), response.position, response.velocity, ...
                response.forceInternalMechanics, response.powerInternalMechanics];
            if any(string(model) == ["RM3_END_STOPS", "RM3_END_STOPS_STEP", ...
                                     "RM3_END_STOPS_STEP_FINE"])
                values = [values, response.forceTotal, ...
                    response.forceConstraint, response.forceActuation, ...
                    response.acceleration];
            end
            assert(all(isfinite(values), 'all'), 'The MATLAB PTO response contains nonfinite values');
            filename = sprintf('%s_%s_pto%d.csv', model, cases(iCase), iPto);
            writematrix(values, fullfile(outDir, filename));
        end
    end
    close_system(erase(simu.simMechanicsFile, '.slx'), 0);
    clear output;
end
end

function [output, simu] = run_radiation_option(option)
% Isolate each published radiation option's script workspace.
enable_convolution = 0;
enable_FIR = 0;
enable_ss = 0;
switch option
    case "constant"
        enable_convolution = 0;
    case "FIR"
        enable_FIR = 1;
    case "state_space"
        enable_ss = 1;
    case "convolution"
        enable_convolution = 1;
    otherwise
        error('Unknown RM3 radiation option: %s', option);
end
wecSim;
end

function run_rm3_mcr_baseline(repoRoot, outDir, model)
% Execute each pinned RM3 MCR driver, including its postprocessing.
switch model
    case "RM3_MCR_ARRAY"
        folder = 'RM3_MCROPT1';
    case "RM3_MCR_EXCEL"
        folder = 'RM3_MCROPT2';
    case "RM3_MCR_MAT"
        folder = 'RM3_MCROPT3';
    otherwise
        error('Unknown RM3 MCR option: %s', model);
end
label = char(model);
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Multiple_Condition_Runs', ...
    folder);
cd(caseDir);
wecSimMCR;
expectedHeader = {'waves.height', 'waves.period', ...
    'pto(1).damping', 'pto(1).stiffness'};
assert(isequal(mcr.header, expectedHeader), ...
    'Pinned RM3 MCR header changed');
assert(isequal(size(mcr.cases), [8, 4]), ...
    'Pinned RM3 MCR must contain eight four-parameter cases');
assert(numel(mcr.Avgpower) == 8 && all(isfinite(mcr.Avgpower)), ...
    'RM3 MCR did not calculate eight finite average powers');
assert(numel(mcr.CPTO) == 8 && all(isfinite(mcr.CPTO)), ...
    'RM3 MCR did not record eight finite PTO damping values');
writematrix([mcr.cases, mcr.Avgpower(:), mcr.CPTO(:)], ...
    fullfile(outDir, sprintf('%s_summary.csv', label)));

% Export the actual CData plotted by the pinned userDefinedFunctionsMCR.m.
% Its axes labels are transposed in the source example, but the tick values
% and 2-by-2 cell placement identify period rows and wave-height columns.
powerImages = findall(groot, 'Type', 'image');
foundD12 = false;
foundD24 = false;
for iImage = 1:numel(powerImages)
    matrix = double(powerImages(iImage).CData);
    if ~isequal(size(matrix), [2, 2])
        continue
    end
    axesHandle = ancestor(powerImages(iImage), 'axes');
    titleText = string(axesHandle.Title.String);
    if contains(titleText, 'Damping = 1200000')
        assert(~foundD12 && all(isfinite(matrix), 'all'), ...
            'RM3 MCR first power matrix is invalid');
        writematrix(matrix, fullfile(outDir, sprintf( ...
            '%s_power_matrix_D12.csv', label)));
        foundD12 = true;
    elseif contains(titleText, 'Damping = 2400000')
        assert(~foundD24 && all(isfinite(matrix), 'all'), ...
            'RM3 MCR second power matrix is invalid');
        writematrix(matrix, fullfile(outDir, sprintf( ...
            '%s_power_matrix_D24.csv', label)));
        foundD24 = true;
    end
end
assert(foundD12 && foundD24, ...
    'RM3 MCR did not plot both published PTO power matrices');

for iCase = 1:8
    saved = load(fullfile(caseDir, sprintf('savedData%03d.mat', iCase)), ...
        'output');
    assert(~isempty(saved.output.bodies) && isstruct(saved.output.ptos), ...
        'RM3 MCR savedData is missing motion or PTO output');
    for iBody = 1:numel(saved.output.bodies)
        response = saved.output.bodies(iBody);
        values = [response.time(:), response.position, response.velocity, ...
            response.forceTotal, response.forceExcitation];
        assert(all(isfinite(values), 'all'), ...
            'RM3 MCR body output contains nonfinite values');
        writematrix(values, fullfile(outDir, sprintf( ...
            '%s_case%d_body%d.csv', label, iCase, iBody)));
    end
    for iPto = 1:numel(saved.output.ptos)
        response = saved.output.ptos(iPto);
        values = [response.time(:), response.position, response.velocity, ...
            response.forceInternalMechanics, response.powerInternalMechanics];
        assert(all(isfinite(values), 'all'), ...
            'RM3 MCR PTO output contains nonfinite values');
        writematrix(values, fullfile(outDir, sprintf( ...
            '%s_case%d_pto%d.csv', label, iCase, iPto)));
    end
end
close_system('RM3', 0);
close all;
end

function run_rm3_mcr_seastate_baseline(repoRoot, outDir)
% Run the three published spectrumImport cases with their imported phases.
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Multiple_Condition_Runs', ...
    'RM3_MCROPT3_SeaState');
cd(caseDir);
wecSimMCR;
assert(isequal(mcr.header, {'waves.spectrumFile', 'simu.solver'}), ...
    'Pinned sea-state MCR header changed');
assert(isequal(size(mcr.cases), [3, 2]), ...
    'Pinned sea-state MCR must contain three cases');
assert(isequal(string(mcr.cases(:,1)), ...
    ["spectrumData1.mat"; "spectrumData2.mat"; "spectrumData3.mat"]), ...
    'Pinned sea-state MCR spectrum order changed');
assert(all(strcmp(mcr.cases(:,2), 'ode4')), ...
    'Pinned sea-state MCR solver changed');
assert(numel(mcr.Avgpower) == 3 && all(isfinite(mcr.Avgpower)), ...
    'Sea-state MCR did not calculate three finite powers');
writematrix(mcr.Avgpower(:), ...
    fullfile(outDir, 'RM3_MCR_SEASTATE_summary.csv'));
% Keep both mass matrices: MATLAB moves part of infinite-frequency added
% mass into the Simscape rigid body before each run, then restores the
% hydrodynamic matrix during postprocessing. The force output alone does
% not reveal which matrix the dynamics used.
for iBody = 1:numel(body)
    hydroForce = body(iBody).hydroForce.hf1;
    assert(isfield(hydroForce.storage, 'hydroForce_fAddedMass'), ...
        'Sea-state MCR is missing the applied added-mass matrix');
    writematrix(hydroForce.fAddedMass, fullfile(outDir, sprintf( ...
        'RM3_MCR_SEASTATE_body%d_added_mass_original.csv', iBody)));
    writematrix(hydroForce.storage.hydroForce_fAddedMass, ...
        fullfile(outDir, sprintf( ...
        'RM3_MCR_SEASTATE_body%d_added_mass_applied.csv', iBody)));
    properties = [body(iBody).mass, hydroForce.mass, ...
        hydroForce.adjustedMass, body(iBody).inertia, ...
        hydroForce.adjustedInertia];
    writematrix(properties, fullfile(outDir, sprintf( ...
        'RM3_MCR_SEASTATE_body%d_mass_properties.csv', iBody)));
end
for iCase = 1:3
    saved = load(fullfile(caseDir, sprintf('savedData%03d.mat', iCase)), ...
        'output', 'waves');
    assert(~isempty(saved.output.bodies) && isstruct(saved.output.ptos), ...
        'Sea-state MCR savedData is missing body or PTO output');
    components = [saved.waves.omega(:), saved.waves.amplitude(:), ...
        saved.waves.dOmega(:), saved.waves.phase(:)];
    assert(size(components, 2) == 4 && all(isfinite(components), 'all'), ...
        'Sea-state MCR has invalid wave components');
    writematrix(components, fullfile(outDir, sprintf( ...
        'RM3_MCR_SEASTATE_case%d_components.csv', iCase)));
    writematrix(saved.waves.waveAmpTime, fullfile(outDir, sprintf( ...
        'RM3_MCR_SEASTATE_case%d_wave.csv', iCase)));
    for iBody = 1:numel(saved.output.bodies)
        response = saved.output.bodies(iBody);
        values = [response.time(:), response.position, response.velocity, ...
            response.forceTotal, response.forceExcitation];
        assert(all(isfinite(values), 'all'), ...
            'Sea-state MCR body output contains nonfinite values');
        writematrix(values, fullfile(outDir, sprintf( ...
            'RM3_MCR_SEASTATE_case%d_body%d.csv', iCase, iBody)));
        forceValues = [response.time(:), response.forceRadiationDamping, ...
            response.forceAddedMass, response.forceRestoring, ...
            response.acceleration];
        assert(all(isfinite(forceValues), 'all'), ...
            'Sea-state MCR hydrodynamic forces contain nonfinite values');
        writematrix(forceValues, fullfile(outDir, sprintf( ...
            'RM3_MCR_SEASTATE_case%d_body%d_forces.csv', iCase, iBody)));
    end
    for iPto = 1:numel(saved.output.ptos)
        response = saved.output.ptos(iPto);
        values = [response.time(:), response.position, response.velocity, ...
            response.forceInternalMechanics, response.powerInternalMechanics];
        assert(all(isfinite(values), 'all'), ...
            'Sea-state MCR PTO output contains nonfinite values');
        writematrix(values, fullfile(outDir, sprintf( ...
            'RM3_MCR_SEASTATE_case%d_pto%d.csv', iCase, iPto)));
    end
end
close_system('RM3', 0);
close all;
end
