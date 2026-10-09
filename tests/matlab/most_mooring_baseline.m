function most_mooring_baseline
% Evaluate the pinned MOST direct nonlinear-static mooring at several poses.
mooring = mooringClass('mooring1');
mooring.nonlinearStaticData = struct(...
    'flag', 1, 'd', 0.333, 'L', 850, 'linearMassAir', 685, ...
    'nLines', 3, 'nodes', [-58 0 -14; -837 0 -inf]', ...
    'EA', 3.27e9, 'CB', 1, 'MaxIter', 150, ...
    'TolFun', 1e-5, 'TolX', 1e-5, ...
    'HV0_try', [1e6 2e6; 1e6 2e6; 1e6 2e6]);
mooring.nonlinearStaticSetup(1025, 9.80665, 200);

degrees = pi/180;
poses = [
    0 0 0 0 0 0;
    5 0 0 0 0 0;
    -5 0 0 0 0 0;
    0 5 0 0 0 0;
    0 0 2.5 0 0 0;
    0 0 -5 0 0 0;
    0 0 0 5*degrees 0 0;
    0 0 0 0 -5*degrees 0;
    0 0 0 0 0 10*degrees;
    10 -5 2.5 5*degrees -5*degrees 10*degrees];
forces = zeros(size(poses,1), 6);
tensions = zeros(size(poses,1), 3, 2);
for i = 1:size(poses,1)
    [f, hv] = nonLinearStaticMooring(poses(i,:), ...
        mooring.nonlinearStaticData.HV0, mooring.nonlinearStaticData);
    forces(i,:) = f';
    tensions(i,:,:) = hv;
end
save('matlab-most-mooring.mat', 'poses', 'forces', 'tensions', '-v7');
end
