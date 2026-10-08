% Run from the Python fork root after adding the pinned MATLAB WEC-Sim paths.
% Produces direct waveClass outputs for comparison with the production port.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = fullfile(repoRoot, 'matlab-reference-output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

times = (0:0.1:200)';
waves = waveClass('regular');
waves.period = 8;
waves.height = 2.5;
waves.direction = 0;
waves.marker.location = [5 5; 10 0; 0 -10];
waves.setup([5.19999512307279, 0.0199999977946844], ...
    'infinite', 100, 0.1, 2000, times, 9.81, 1000);
writematrix(waves.waveAmpTime, fullfile(outputDir, 'regular_origin.csv'));
writematrix(waves.waveAmpTimeViz, fullfile(outputDir, 'regular_markers.csv'));

finite = waveClass('regular');
finite.period = 8;
finite.height = 2.5;
finite.direction = 0;
finite.setup([0.4, 2.0], 25, 0, 0.1, 2000, times, 9.81, 1000);
writematrix([finite.wavenumber, finite.power], ...
    fullfile(outputDir, 'finite_depth.csv'));

% Current IEC PM and JONSWAP spectra. Only deterministic spectral quantities
% are paired here: MATLAB's seeded phases use Threefry substreams.
for spectrum = ["PM", "JS"]
    irregular = waveClass('irregular');
    irregular.period = 8;
    irregular.height = 2.5;
    if spectrum == "JS"
        irregular.height = 4;
    end
    irregular.spectrumType = char(spectrum);
    irregular.bem.option = 'Traditional';
    irregular.bem.count = 64;
    irregular.phaseSeed = 1;
    irregular.setup([0.4, 2.0], 'infinite', 1, 0.1, 20, ...
        (0:0.1:2)', 9.81, 1000);
    label = lower(char(spectrum));
    writematrix([irregular.omega, irregular.spectrum], ...
        fullfile(outputDir, sprintf('%s_spectrum.csv', label)));
    writematrix(irregular.power, ...
        fullfile(outputDir, sprintf('%s_power.csv', label)));
    if spectrum == "JS"
        writematrix(irregular.gamma, ...
            fullfile(outputDir, 'js_gamma.csv'));
    end
end

fid = fopen(fullfile(outputDir, 'matlab_release.txt'), 'w');
fprintf(fid, '%s\n', version('-release'));
fclose(fid);
