function reference_model_baseline(model)
% Run the published RM3, OSWEC, or Sphere reference cases and save numeric output.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
outDir = fullfile(repoRoot, 'matlab-reference-model-output');
if ~isfolder(outDir)
    mkdir(outDir);
end

switch string(model)
    case {"RM3", "OSWEC"}
        cases = string(model);
        caseDirs = string(fullfile(repoRoot, 'matlab-ref', 'examples', model));
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
        cases = ["B2B_Case1", "B2B_Case2"];
        caseDirs = fullfile(repoRoot, 'applications', ...
            'Body-to-Body_Interactions', cases);
    otherwise
        error('Unknown reference model: %s', model);
end

for iCase = 1:numel(cases)
    cd(caseDirs(iCase));
    wecSim;
    assert(exist('output', 'var') == 1 && ~isempty(output.bodies), ...
        'The MATLAB case produced no body output');
    if string(model) == "OSWEC"
        % The published case shuffles its random phase on each run. Save the
        % actual components so wave and force checks use the paired realization.
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
        assert(all(isfinite(values), 'all'), 'The MATLAB response contains nonfinite values');
        filename = sprintf('%s_%s_body%d.csv', model, cases(iCase), iBody);
        writematrix(values, fullfile(outDir, filename));
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
