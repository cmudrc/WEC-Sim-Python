function application_regression(folder)
% Execute the upstream WEC-Sim Applications tests for a reference-model folder.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
appRoot = fullfile(repoRoot, 'applications');
cd(appRoot);
if ismember(string(folder), ["Mooring", "Paraview_Visualization"])
    % Upstream skips MoorDyn cases solely when GITHUB_ACTIONS is true.
    originalCI = getenv('GITHUB_ACTIONS');
    restoreCI = onCleanup(@() setenv('GITHUB_ACTIONS', originalCI)); %#ok<NASGU>
    setenv('GITHUB_ACTIONS', 'false');
end
results = wecSimAppTest(string(folder));
outDir = fullfile(repoRoot, 'matlab-application-results');
if ~isfolder(outDir)
    mkdir(outDir);
end
if isempty(results)
    caseDir = fullfile(appRoot, folder);
    assert(isfile(fullfile(caseDir, 'wecSimInputFile.m')), ...
        'No upstream tests or runnable input file for %s', folder);
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
