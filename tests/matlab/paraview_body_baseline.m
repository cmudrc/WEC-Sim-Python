function paraview_body_baseline
% Exercise the pinned body VTP writer with translation, XYZ rotation, and pressures.
body.geometry.vertex = [0 0 0; 2 0 0; 0 1 0; 0 0 2];
body.geometry.face = [1 2 3; 1 3 4];
body.geometry.area = [1; 1];
body.geometry.numVertex = 4;
body.geometry.numFace = 2;
times = [2; 3];
poses = [0 0 0 0 0 0; 1 -2 3 0.2 -0.3 0.4];
hydrostatic = [10 20; 30 40];
nonlinear = [1 2; 3 4];
linear = [5 6; 7 8];
output = fullfile(pwd, 'matlab-paraview-body');
mkdir(output);
mkdir(fullfile(output, 'body1_flap'));
writeParaviewBody(body, times, poses, 'flap', 'OSWEC', 'pinned', ...
    hydrostatic, nonlinear, linear, output, 1);
end
