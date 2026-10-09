function paraview_wave_baseline
% Export two small regular-wave grids with the pinned WEC-Sim VTP writer.
waves = waveClass('regular');
waves.height = 2.5;
waves.period = 8;
waves.direction = 90;
times = [2 3];
waves.setup([0.3 1.5], 30, 0, 1, length(times), times, 9.81, 1000);
output = fullfile(pwd, 'matlab-paraview-wave');
mkdir(output);
mkdir(fullfile(output, 'waves'));
writeParaviewWave(waves, times, 3, 2, 6, 'OSWEC', 'pinned', ...
    0, output, times, 9.81);
writematrix([waves.amplitude waves.wavenumber waves.omega ...
    waves.direction waves.waterDepth], fullfile(output, 'parameters.csv'));
end
