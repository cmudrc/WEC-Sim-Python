function application_regression(folder)
% Execute the upstream WEC-Sim Applications tests for a reference-model folder.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
cd(fullfile(repoRoot, 'applications'));
results = wecSimAppTest(string(folder));
outDir = fullfile(repoRoot, 'matlab-application-results');
if ~isfolder(outDir)
    mkdir(outDir);
end
summary = table(string({results.Name})', [results.Passed]', ...
    [results.Failed]', [results.Incomplete]', [results.Duration]', ...
    'VariableNames', {'Name', 'Passed', 'Failed', 'Incomplete', 'Duration'});
writetable(summary, fullfile(outDir, sprintf('%s.csv', folder)));
assertSuccess(results);
end
