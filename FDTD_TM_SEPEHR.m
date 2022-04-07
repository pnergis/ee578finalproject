%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sepehr Eskandari
% EE578: Computational Electromagnetics for Engineers, Spring 2022
% Date: 3/16/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars;
clc;

% Domain Size
xDim = 5000e-9;
% Grid step size
dx = 20e-9;
yDim = 5000e-9;
dy = 20e-9;
% Time step
dt = dx/sqrt(2);
% Total number of time steps
nTotal = 10e3;
% Number of grids in each direction
Nx = floor(xDim/dx);
Ny = floor(yDim/dy);
% Inhomogeneity
er = ones(Ny+1,Nx+1);
er(66:187,66:187) = 12;
er(151:173,186:end) = 12;
er(79:101,186:end) = 12;
er(115:137,1:65) = 12;
% Initialization of fields
Ez = zeros(Ny+1,Nx+1);
Dz = zeros(Ny+1,Nx+1);
Hx = zeros(Ny,Nx+1);
Bx = zeros(Ny,Nx+1);
Hy = zeros(Ny+1,Nx);
By = zeros(Ny+1,Nx);
% Source parameters
lambdaL = 1534e-9;
lambda0 = 1550e-9;
lambdaU = 1566e-9;
w0 = 2*pi/lambda0;
sigma = 2*lambda0/w0/(lambdaU-lambdaL);
% Source and probe locations
%xSource = -1850e-9;
%ySource = 0;
x_observation = 2020e-9;
y_observation = 500e-9;
[X,Y] = meshgrid(linspace(-xDim/2,xDim/2,Nx+1),linspace(-yDim/2,yDim/2,Ny+1));
%diff_source = sqrt((X-xSource).^2+(Y-ySource).^2);
%[~,idx]=min(diff_source(:));
%[row_source,col_source]=ind2sub(size(diff_source),idx);
diff_observation = sqrt((X-x_observation).^2+(Y-y_observation).^2);
[~,idx]=min(diff_observation(:));
[row_observation,col_observation]=ind2sub(size(diff_observation),idx);
% PML width cell
L = 15;
% Gradung order
gradingOrder = 3;
% Maximum reflection
maxReflection = 1e-9;
% Maximum sigma
maxSigma = -(gradingOrder+1)/4*(1/L/dx)*log(maxReflection);
% SigmaY evaluated at Bx location
% Bx1 and Bx2 are multiplication factor matrices.
Bx_SigmaY = zeros(Ny,Nx+1);
for i=1:L
    Bx_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L).^gradingOrder;
    Bx_SigmaY(Ny-i+1,:) = maxSigma*(1+(1/2-i)/L).^gradingOrder;
end
Bx1 = (1-Bx_SigmaY*dt/2)./(1+Bx_SigmaY*dt/2);
Bx2 = 1+Bx_SigmaY*dt/2;


% SigmaX evaluated at By location
% By1 and By2 are multiplication factor matrices.
By_SigmaX = zeros(Ny+1,Nx);
for i=1:L
    By_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    By_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
By1 = (1-By_SigmaX*dt/2)./(1+By_SigmaX*dt/2);
By2 = 1+By_SigmaX*dt/2;


% SigmaX evaluated at Hx location
% Hx1 and Hx2 are multiplication factor matrices.
Hx_SigmaX = zeros(Ny,Nx+1);
for i=1:L
    Hx_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Hx_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Hx1 = 1 + Hx_SigmaX*dt/2;
Hx2 = 1 - Hx_SigmaX*dt/2;


% SigmaY evaluated at Hy location
% Hy1 and Hy2 are multiplication factor matrices.
Hy_SigmaY = zeros(Ny+1,Nx);
for i=1:L
    Hy_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Hy_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Hy1 = 1 + Hy_SigmaY*dt/2;
Hy2 = 1 - Hy_SigmaY*dt/2;


% SigmaX evaluated at Dz location
% Dz1 and Dz2 are multiplication factor matrices.
Dz_SigmaX = zeros(Ny+1,Nx+1);
for i=1:L
    Dz_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Dz_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Dz1 = (1-Dz_SigmaX*dt/2)./(1+Dz_SigmaX*dt/2);
Dz2 = 1+Dz_SigmaX*dt/2;


% SigmaY evaluated at Ez location.
% Ez1 and Ez2 are multiplication factor matrices.
Ez_SigmaY = zeros(Ny+1,Nx+1);
for i=1:L
   Ez_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
   Ez_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Ez1 = (1-Ez_SigmaY*dt/2)./(1+Ez_SigmaY*dt/2);
Ez2 = 1+Ez_SigmaY*dt/2;
% For defining a figure of merit, I recored electric field at
% (2020nm,500nm) as the observation point and at (0,480nm) as excitation 
% point and then calculated the ratio of maximum of probe location field to
% the excitation location field.
Ez_Probe = zeros(nTotal,1);
Ez_Excitation = zeros(nTotal,1);
Ez_2DFigure = zeros(Ny+1,Nx+1);
% Update loop
%J = zeros(nTotal,1);
for n=1:nTotal
   Bx_BeforeUpdate = Bx;
   By_BeforeUpdate = By;
   Dz_BeforeUpdate = Dz;
   Bx(:,2:Nx) = Bx1(:,2:Nx).*Bx(:,2:Nx) - dt./dy./Bx2(:,2:Nx).*(Ez(2:end,2:Nx)-Ez(1:Ny,2:Nx)); % step 1
   By(2:Ny,:) = By1(2:Ny,:).*By(2:Ny,:) + dt./dx./By2(2:Ny,:).*(Ez(2:Ny,2:end)-Ez(2:Ny,1:Nx)); % step 2
   Hx = Hx + Hx1.*Bx - Hx2.*Bx_BeforeUpdate; % step 3
   Hy = Hy + Hy1.*By - Hy2.*By_BeforeUpdate; % step 4
   Dz(2:Ny,2:Nx) = Dz1(2:Ny,2:Nx).*Dz(2:Ny,2:Nx) + dt./Dz2(2:Ny,2:Nx)./er(2:Ny,2:Nx).*((Hy(2:Ny,2:end)-Hy(2:Ny,1:(Nx-1)))/dx-(Hx(2:end,2:Nx)-Hx(1:(Ny-1),2:Nx))/dy); % step 5
   Ez(2:Ny,2:Nx) = Ez1(2:Ny,2:Nx).*Ez(2:Ny,2:Nx) + 1./Ez2(2:Ny,2:Nx).*(Dz(2:Ny,2:Nx)-Dz_BeforeUpdate(2:Ny,2:Nx)); % step 6
   %%%%% Fundamental Odd TM Mode
%    Ez(127:137,L+10) = exp(-0.8936*2.2)/sin(1.0308*2.2)*sin(1.0308e7*(1:11)*dy)*cos(w0*n*dt);
%    Ez(126,L+10) = 0;
%    Ez(115:125,L+10) = -flipud(Ez(127:137,L+10));
%    Ez(138:220,L+10) = exp(-0.8936*1e7*dy*(12:94))*cos(w0*n*dt);
   Ez(32:114,L+10) = -flipud(Ez(138:220,L+10));
   %%%%% Fundamental Even TM Mode
   % Right Now the sources are hard. Need to be modified later.
   Ez(127:137,L+10) = (exp(-1.236*2.2)/cos(0.5299*2.2)*cos(5.299e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
   Ez(126,L+10) = (exp(-1.236*2.2)/cos(0.5299*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
   Ez(115:125,L+10) = flipud(Ez(127:137,L+10));
   Ez(138:220,L+10) = (exp(-1.236*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
   Ez(32:114,L+10) = flipud(Ez(138:220,L+10));

%    %%%%% Fundamental Even TM Mode
%    Ez(127:137,L+10) = Ez(127:137,L+10) + (exp(-1.236*2.2)/cos(0.5299*2.2)*cos(5.299e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
%    Ez(126,L+10) = Ez(126,L+10) + (exp(-1.236*2.2)/cos(0.5299*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
%    Ez(115:125,L+10) = flipud(Ez(127:137,L+10));
%    Ez(138:220,L+10) = Ez(138:220,L+10) + (exp(-1.236*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
%    Ez(32:114,L+10) = flipud(Ez(138:220,L+10));
%     

   %%%%%
   Ez_Probe(n) = Ez(row_observation,col_observation);
   Ez_Excitation(n) = Ez(126,L+10);
%    if(n>4000)
%        pcolor(X,Y,Ez);
%        title(['Electric Field Profile at ', num2str(n), 'th time-step'])
%        colorbar
%        rectangle('Position',[-1200e-9,-1200e-9,2400e-9,2400e-9]);
%        rectangle('Position',[-2500e-9,-220e-9,1300e-9,440e-9]);
%        rectangle('Position',[1200e-9,-720e-9,1300e-9,440e-9]);
%        rectangle('Position',[1200e-9,280e-9,1300e-9,440e-9]);
%        caxis([-0.1 0.1])
%        shading interp
%        pause(0.01);
%    end
n
end
efficiency = max(Ez_Probe) / max(Ez_Excitation);