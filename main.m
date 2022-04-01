% % EE 578 Final Project
% Sepehr Eskendari and Piril Nergis
% main.m calls the optimizer calls FDTD_2D_TM.m and FDTD_2D_TE as it
% iterates/optimizes
%
% The polarization beamsplitter is 2.4 micron by 2.4 micron and is
% represented by 20x20 pixels

% Polarization Beamsplitter (pbs)

%generates initial 20x20 matrix for optimization
Npix=20;
pbs = randi([0 1], Npix,Npix);

[Ez, Hx, Hy] = FDTD_2D_TM(pbs);