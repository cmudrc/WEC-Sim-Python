function threefry_reference()
% Inspect the pinned WEC-Sim waveClass phase stream without a WEC model.
outDir = fullfile(pwd, 'matlab-threefry-output');
if ~isfolder(outDir)
    mkdir(outDir);
end
for substream = 1:3
    stream = RandStream('Threefry', 'Seed', 1);
    stream.Substream = substream;
    before = stream.State;
    samples = rand(stream, 1, 32);
    after = stream.State;
    save(fullfile(outDir, sprintf('seed1_substream%d.mat', substream)), ...
        'before', 'after', 'samples', 'substream');
    writematrix(samples, fullfile(outDir, ...
        sprintf('seed1_substream%d.csv', substream)));
    fprintf('Threefry seed 1 substream %d: state %s [%s], first rand %.17g\n', ...
        substream, class(before), num2str(size(before)), samples(1));
end
end
