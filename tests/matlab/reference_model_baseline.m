function reference_model_baseline(model)
% Run the published RM3, OSWEC, or Sphere reference cases and save numeric output.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
outDir = fullfile(repoRoot, 'matlab-reference-model-output');
if ~isfolder(outDir)
    mkdir(outDir);
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
        if string(model) == "RM3_Radiation_Options"
            values = [values, response.forceRadiationDamping, ...
                response.forceAddedMass];
        end
        assert(all(isfinite(values), 'all'), 'The MATLAB response contains nonfinite values');
        filename = sprintf('%s_%s_body%d.csv', model, cases(iCase), iBody);
        writematrix(values, fullfile(outDir, filename));
    end
    if ismember(string(model), ["Sphere_Passive", "Sphere_PTO_Config", ...
            "Sphere_Reactive_PI"])
        assert(exist('controller1_out', 'var') == 1, ...
            'The passive controller produced no logged output');
        controllerValues = [controller1_out.time(:), controller1_out.signals.values];
        assert(size(controllerValues, 2) == 13, ...
            'Expected six force and six power components from passive controller');
        assert(all(isfinite(controllerValues), 'all'), ...
            'The passive controller output contains nonfinite values');
        writematrix(controllerValues, fullfile(outDir, string(model) + "_controller.csv"));
    end
    if isstruct(output.ptos) && isfield(output.ptos, 'time')
        for iPto = 1:numel(output.ptos)
            response = output.ptos(iPto);
            values = [response.time(:), response.position, response.velocity, ...
                response.forceInternalMechanics, response.powerInternalMechanics];
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
