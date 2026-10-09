function most_bem_baseline
% Evaluate the pinned MOST library BEM block on still and moving hub states.
root = pwd;
source = fullfile(root, 'applications', 'MOST', 'mostData', ...
    'windTurbine', 'turbine_properties');
model = windTurbineClass('IEA15MW');
model.turbineName = fullfile(source, 'Properties_IEA15MW.mat');
model.bladeDataName = fullfile(source, 'Bladedata_IEA15MW.mat');
model.loadTurbineData();
x = -30:10:60;
y = linspace(-150, 150, 12);
z = linspace(10, 290, 12);
model.createBEMstruct(x, y, z);
bem = model.BEMstruct;

q = zeros(4, 14);
q(:,3) = model.hub.height;
q(:,14) = 7.56*pi/30;
q(2,13) = pi/6;
q(3,1:3) = [5 2 model.hub.height+1];
q(3,4:6) = [3 -2 5]*pi/180;
q(3,7:9) = [0.1 -0.05 0.03];
q(3,10:12) = [0.002 -0.003 0.001];
q(3,13) = pi/4;
q(4,1:3) = [-3 -2 model.hub.height-2];
q(4,4:6) = [-2 3 -4]*pi/180;
q(4,7:9) = [-0.1 0.08 -0.04];
q(4,10:12) = [-0.002 0.001 -0.003];
q(4,13) = pi/2;
q(4,14) = 0.9;
wind_speeds = [8 8 8 11];
bladepitch = [0 0 0.03 0.08];
loads = zeros(4, 6, 3);
wind = zeros(length(x), length(y), length(z), 3);
for i = 1:4
    wind(:,:,:,1) = wind_speeds(i);
    loads(i,:,:) = BEM(wind, q(i,:)', bladepitch(i), bem);
end
assert(all(isfinite(loads), 'all'));
save('matlab-most-bem.mat', 'bem', 'q', 'wind_speeds', ...
    'bladepitch', 'loads', '-v7');
end
