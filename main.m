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
writematrix(['FOM_TE' num2str(FOM_TE); 'FOM_TM' num2str(FOM_TM), '/project/mahta_676/nergis/ee578finalproject/results/04_14_1pm_FOM_results.txt'])
writematrix(pbs, '/project/mahta_676/nergis/ee578finalproject/results/04_14_1pm_finalPBS.txt','Delimiter','tab')
