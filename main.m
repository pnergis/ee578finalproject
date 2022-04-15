% % EE 578 Final Project
% Sepehr Eskendari and Piril Nergis
% main.m calls the optimizer calls FDTD_2D_TM.m and FDTD_2D_TE as it
% iterates/optimizes
%
% The polarization beamsplitter is 2.4 micron by 2.4 micron and is
% represented by 20x20 pixels

% Polarization Beamsplitter (pbs)
% generates initial 20x20 matrix for optimization
Npix=20;
% pbs = randi([0 1], Npix,Npix);
pbs = ones(20,20);
tic
[FOM_TE, FOM_TM, erTE, erTM, pbs] = optimizeDBS(pbs);
toc

%% validate optimized pbs
[efficiency, ~, Ez, Hy, Hx] = FDTD_2D_TM(pbs);


%%
figure (1)
pcolor(-2500e-9:20e-9:2500e-9,-2500e-9:20e-9:2500e-9,erTM)
shading interp
xlabel('x')
ylabel('y')
title('\epsilon of Domain');
set(gca,'YDir','normal')

%%
pbs = randi([0 1], Npix,Npix);
    % Domain Size
    xDim = 5000e-9;
    yDim = 5000e-9;
    % Grid step size
    dx = 20e-9;
    dy = 20e-9;
    % Time step
    dt = dx/sqrt(2);
    % Total number of time steps
    nTotal = 10e3;
    % Number of grids in each direction
    Nx = floor(xDim/dx);
    Ny = floor(yDim/dy);
    % PBS pixel to epsilon grid
    er = ones(Ny+1,Nx+1);
    erSi = 12;
    wg_size = 440e-9;
    wg_dist = 1000e-9; % distance between the center of output wgs
    pbs_size = 2400e-9;
%   er(151:173,186:end) = 12;
% 	er(79:101,186:end) = 12;
    er((Ny/2-wg_size/dy/2)+1:(Ny/2+wg_size/dy/2)+1,...
        1:(Nx/2-pbs_size/dx/2)) = erSi;
    er(((Ny/2+wg_dist/dy/2)-wg_size/dy/2):((Ny/2+wg_dist/dy/2)+wg_size/dy/2),...
        (Nx/2+pbs_size/dx/2)+1:end) = erSi;
    er(((Ny/2-wg_dist/dy/2)-wg_size/dy/2):((Ny/2-wg_dist/dy/2)+wg_size/dy/2),...
        (Nx/2+pbs_size/dx/2)+1:end) = erSi;
    er = pixelToEr(pbs, er, erSi, Ny, Nx, dy, dx, pbs_size); % pbs
figure (1)
pcolor(-2500e-9:20e-9:2500e-9,-2500e-9:20e-9:2500e-9,er)
xlabel('x')
ylabel('y') 
title('\epsilon of Domain');
set(gca,'YDir','normal')