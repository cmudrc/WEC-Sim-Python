function most_bem_advection_baseline
% Feed published TurbSim time-shifted wind frames into the active BEM block.
rotor = load('matlab-most-bem.mat', 'bem', 'q', 'bladepitch');
advection = load('matlab-most-advection.mat', 'values', 'x', 'y', 'z');
bem = rotor.bem;
bem.WIND_Xdiscr = advection.x;
bem.WIND_Ydiscr = advection.y;
bem.WIND_Zdiscr = advection.z;
indices = [1 125 249];
q = rotor.q([1 3 4], :);
bladepitch = rotor.bladepitch([1 3 4]);
points = [0 0 150; 5 35 220; -3 -50 80];
loads = zeros(3, 6, 3);
sampled_wind = zeros(3, 3, 3);
for case_index = 1:3
    frame = squeeze(advection.values(indices(case_index), :, :, :, :));
    wind = permute(frame, [2 3 4 1]);
    loads(case_index, :, :) = BEM(wind, q(case_index, :)', ...
        bladepitch(case_index), bem);
    for point_index = 1:3
        point = points(point_index, :);
        for component = 1:3
            sampled_wind(case_index, point_index, component) = interp3( ...
                advection.y, advection.x, advection.z, wind(:, :, :, component), ...
                point(2), point(1), point(3));
        end
    end
end
assert(all(isfinite(loads), 'all'));
assert(all(isfinite(sampled_wind), 'all'));
save('matlab-most-bem-advection.mat', 'indices', 'q', 'bladepitch', ...
    'points', 'sampled_wind', 'loads', '-v7');
end
