function application_regression(folder, forceMoorDyn)
% Execute the upstream WEC-Sim Applications tests for a reference-model folder.
if nargin < 2
    forceMoorDyn = 'false';
end
if strcmpi(forceMoorDyn, 'true') && ...
        any(string(folder) == ["Mooring", "Paraview_Visualization"])
    originalActions = getenv('GITHUB_ACTIONS');
    restoreActions = onCleanup(@() setenv('GITHUB_ACTIONS', originalActions));
    setenv('GITHUB_ACTIONS', 'false');
end
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
appRoot = fullfile(repoRoot, 'applications');
cd(appRoot);
results = wecSimAppTest(string(folder));
outDir = fullfile(repoRoot, 'matlab-application-results');
if ~isfolder(outDir)
    mkdir(outDir);
end
if isempty(results)
    caseDir = fullfile(appRoot, folder);
    assert(isfile(fullfile(caseDir, 'wecSimInputFile.m')), ...
        'No upstream tests or runnable input file for %s', folder);
    if string(folder) == "Multiple_Wave_Spectra"
        hydroDir = fullfile(appRoot, '_Common_Input_Files', 'OSWEC', 'hydroData');
        cd(hydroDir);
        if ~isfile('oswec.h5')
            bemio;
        end
    end
    cd(caseDir);
    wecSim;
    assert(exist('output', 'var') == 1 && ~isempty(output.bodies), ...
        'The direct application case produced no body output');
    summary = table("direct/wecSim", true, false, false, NaN, ...
        'VariableNames', {'Name', 'Passed', 'Failed', 'Incomplete', 'Duration'});
else
    summary = table(string({results.Name})', [results.Passed]', ...
        [results.Failed]', [results.Incomplete]', [results.Duration]', ...
        'VariableNames', {'Name', 'Passed', 'Failed', 'Incomplete', 'Duration'});
end
writetable(summary, fullfile(outDir, sprintf('%s.csv', folder)));
if ~isempty(results)
    assertSuccess(results);
end
end
