function wavestar_nmpc_forecast_probe
% Isolate the published AR/forecast call just after NMPC activates.
outDir = fullfile(pwd, 'matlab-reference-model-output');
saved = load(fullfile(outDir, ...
    'WECCCOMP_NMPC_FINE_ACTIVE_DIAG_controller.mat'));
time = saved.estimated_states.time(:);
estimate = saved.estimated_states.signals.values(:,5);
logged = reshape(saved.AR_excM_pred.signals.values, 40, []);
assert(numel(time) == 15201 && numel(estimate) == 15201 && ...
    isequal(size(logged), [40, 15201]), ...
    'The fine-step NMPC forecast trace changed');

fitIndex = 15001;  % t = 15 s, where the pinned predictor refits.
model = ar(estimate(fitIndex-179:fitIndex), 18);
indices = [15001, 15002, 15011, 15051, 15101, 15201];
rows = zeros(numel(indices), 81);
for j = 1:numel(indices)
    index = indices(j);
    past = estimate(index-179:index);
    direct = forecast(model, past, 40);
    rows(j,:) = [time(index), direct(:).', logged(:,index).'];
end
writematrix(rows, fullfile(outDir, ...
    'WECCCOMP_NMPC_FINE_ACTIVE_DIAG_forecast_probe.csv'));
[A, ~, ~, ~, ~] = polydata(model);
writematrix(A, fullfile(outDir, ...
    'WECCCOMP_NMPC_FINE_ACTIVE_DIAG_ar_coefficients.csv'));
end
