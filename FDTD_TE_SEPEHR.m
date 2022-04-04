%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sepehr Eskandari
% EE578: Computational Electromagnetics for Engineers, Spring 2022
% Date: 3/16/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
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
nTotal = 1e4;
% Number of grids in each direction
Nx = floor(xDim/dx);
Ny = floor(yDim/dy);
% Inhomogeneity

er = ones(Ny+1,Nx+1);
er_Dx = ones(Ny+1,Nx);
er_Dy = ones(Ny,Nx+1);
% er(151:173,186:end) = 12;
% er(79:101,186:end) = 12;
% er(115:137,1:65) = 12;
% er(66:186,66:186) = double(12*(rand(121,121) < 0.5));
% er(66:186,66:186) = 12;
% Initialization of fields
Ez = zeros(Ny+1,Nx+1);
Dz = zeros(Ny+1,Nx+1);
Ex = zeros(Ny+1,Nx);
Dx = zeros(Ny+1,Nx);
Ey = zeros(Ny,Nx+1);
Dy = zeros(Ny,Nx+1);
Hx = zeros(Ny,Nx+1);
Bx = zeros(Ny,Nx+1);
Hy = zeros(Ny+1,Nx);
By = zeros(Ny+1,Nx);
Hz = zeros(Ny,Nx);
Bz = zeros(Ny,Nx);
% Source parameters
lambdaL = 1534e-9;
lambda0 = 1550e-9;
lambdaU = 1566e-9;
w0 = 2*pi/lambda0;
sigma = 2*lambda0/w0/(lambdaU-lambdaL);
% Source and probe locations
xSource = -1850e-9;
ySource = 0;
x_observation = 1000e-9;
y_observation = 0;
[X,Y] = meshgrid(linspace(-xDim/2,xDim/2,Nx+1),linspace(-yDim/2,yDim/2,Ny+1));
diff_source = sqrt((X-xSource).^2+(Y-ySource).^2);
[~,idx]=min(diff_source(:));
[row_source,col_source]=ind2sub(size(diff_source),idx);
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

% SigmaY evaluated at Dx location
% Dx1 and Dx2 are multiplication factor matrices.
Dx_SigmaY = zeros(Ny+1,Nx);
for i=1:L
    Dx_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Dx_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Dx1 = (1-Dx_SigmaY*dt/2)./(1+Dx_SigmaY*dt/2);
Dx2 = 1+Dx_SigmaY*dt/2;


% SigmaX evaluated at By location
% By1 and By2 are multiplication factor matrices.
By_SigmaX = zeros(Ny+1,Nx);
for i=1:L
    By_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    By_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
By1 = (1-By_SigmaX*dt/2)./(1+By_SigmaX*dt/2);
By2 = 1+By_SigmaX*dt/2;


% SigmaX evaluated at By location
% Dy1 and Dy2 are multiplication factor matrices.
Dy_SigmaX = zeros(Ny,Nx+1);
for i=1:L
    Dy_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Dy_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Dy1 = (1-Dy_SigmaX*dt/2)./(1+Dy_SigmaX*dt/2);
Dy2 = 1+Dy_SigmaX*dt/2;


% SigmaX evaluated at Hx location
% Hx1 and Hx2 are multiplication factor matrices.
Hx_SigmaX = zeros(Ny,Nx+1);
for i=1:L
    Hx_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Hx_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Hx1 = 1 + Hx_SigmaX*dt/2;
Hx2 = 1 - Hx_SigmaX*dt/2;


% SigmaX evaluated at Ex location
% Ex1 and Ex2 are multiplication factor matrices.
Ex_SigmaX = zeros(Ny+1,Nx);
for i=1:L
    Ex_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Ex_SigmaX(:,Nx+1-i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Ex1 = 1 + Ex_SigmaX*dt/2;
Ex2 = 1 - Ex_SigmaX*dt/2;


% SigmaY evaluated at Hy location
% Hy1 and Hy2 are multiplication factor matrices.
Hy_SigmaY = zeros(Ny+1,Nx);
for i=1:L
    Hy_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Hy_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Hy1 = 1 + Hy_SigmaY*dt/2;
Hy2 = 1 - Hy_SigmaY*dt/2;


% SigmaY evaluated at Ey location
% Ey1 and Ey2 are multiplication factor matrices.
Ey_SigmaY = zeros(Ny,Nx+1);
for i=1:L
    Ey_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Ey_SigmaY(Ny+1-i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Ey1 = 1 + Ey_SigmaY*dt/2;
Ey2 = 1 - Ey_SigmaY*dt/2;


% SigmaX evaluated at Dz location
% Dz1 and Dz2 are multiplication factor matrices.
Dz_SigmaX = zeros(Ny+1,Nx+1);
for i=1:L
    Dz_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Dz_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Dz1 = (1-Dz_SigmaX*dt/2)./(1+Dz_SigmaX*dt/2);
Dz2 = 1+Dz_SigmaX*dt/2;


% SigmaX evaluated at Bz location
% Bz1 and Bz2 are multiplication factor matrices.
Bz_SigmaX = zeros(Ny,Nx);
for i=1:L
    Bz_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Bz_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Bz1 = (1-Bz_SigmaX*dt/2)./(1+Bz_SigmaX*dt/2);
Bz2 = 1+Bz_SigmaX*dt/2;


% SigmaY evaluated at Ez location.
% Ez1 and Ez2 are multiplication factor matrices.
Ez_SigmaY = zeros(Ny+1,Nx+1);
for i=1:L
   Ez_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
   Ez_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Ez1 = (1-Ez_SigmaY*dt/2)./(1+Ez_SigmaY*dt/2);
Ez2 = 1+Ez_SigmaY*dt/2;


% SigmaY evaluated at Hz location.
% Hz1 and Hz2 are multiplication factor matrices.
Hz_SigmaY = zeros(Ny,Nx);
for i=1:L
   Hz_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
   Hz_SigmaY(Ny+1-i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Hz1 = (1-Hz_SigmaY*dt/2)./(1+Hz_SigmaY*dt/2);
Hz2 = 1+Hz_SigmaY*dt/2;


Ez_Probe = zeros(nTotal,1);
Ez_2DFigure = zeros(Ny+1,Nx+1);
% Update loop
J = zeros(nTotal,1);
for n=1:nTotal
   Bx_BeforeUpdate = Bx;
   Dx_BeforeUpdate = Dx;
   By_BeforeUpdate = By;
   Dy_BeforeUpdate = Dy;
   Dz_BeforeUpdate = Dz;
   Bz_BeforeUpdate = Bz;
   Bx(:,2:Nx) = Bx1(:,2:Nx).*Bx(:,2:Nx) - dt./dy./Bx2(:,2:Nx).*(Ez(2:end,2:Nx)-Ez(1:Ny,2:Nx)); % step 1 - TM
   Dx(2:Ny,:) = Dx1(2:Ny,:).*Dx(2:Ny,:) + dt./dy./Dx2(2:Ny,:)./er_Dx(2:Ny,:).*(Hz(2:Ny,:)-Hz(1:(Ny-1),:)); % step 1 - TE
   By(2:Ny,:) = By1(2:Ny,:).*By(2:Ny,:) + dt./dx./By2(2:Ny,:).*(Ez(2:Ny,2:end)-Ez(2:Ny,1:Nx)); % step 2 - TM
   Dy(:,2:Nx) = Dy1(:,2:Nx).*Dy(:,2:Nx) - dt./dx./Dy2(:,2:Nx)./er_Dy(:,2:Nx).*(Hz(:,2:Nx)-Hz(:,1:(Nx-1))); % step 2 - TE
   Hx = Hx + Hx1.*Bx - Hx2.*Bx_BeforeUpdate; % step 3 - TM
   Ex = Ex + Ex1.*Dx - Ex2.*Dx_BeforeUpdate; % step 3 - TE
   Hy = Hy + Hy1.*By - Hy2.*By_BeforeUpdate; % step 4 - TM
   Ey = Ey + Ey1.*Dy - Ey2.*Dy_BeforeUpdate; % step 4 - TE 
   Dz(2:Ny,2:Nx) = Dz1(2:Ny,2:Nx).*Dz(2:Ny,2:Nx) + dt./Dz2(2:Ny,2:Nx)./er(2:Ny,2:Nx).*((Hy(2:Ny,2:end)-Hy(2:Ny,1:(Nx-1)))/dx-(Hx(2:end,2:Nx)-Hx(1:(Ny-1),2:Nx))/dy); % step 5 - TM
   Bz = Bz1.*Bz + dt./Bz2.*((Ex(2:end,:)-Ex(1:Ny,:))./dy-(Ey(:,2:end)-Ey(:,1:Nx))./dx); % step 5 - TE
   Ez(2:Ny,2:Nx) = Ez1(2:Ny,2:Nx).*Ez(2:Ny,2:Nx) + 1./Ez2(2:Ny,2:Nx).*(Dz(2:Ny,2:Nx)-Dz_BeforeUpdate(2:Ny,2:Nx)); % step 6 - TM
   Hz = Hz1.*Hz + 1./Hz2.*(Bz-Bz_BeforeUpdate);
%    J(n) = exp(-((n-1/2)*dt-4*sigma)^2/(sigma)^2)*sin(w0*((n-1/2)*dt-4*sigma));
   J(n) = sin(w0*((n-1/2)*dt-4*sigma));
   Ez(row_source,col_source) = Ez(row_source,col_source) - dt/er(row_source,col_source)*J(n);
   Hz(row_source,col_source) = Hz(row_source,col_source) - dt*J(n);
   Ez_Probe(n) = Ez(row_observation,col_observation);
%    if(n>3000)
   pcolor(Hz);
   title(['Electric Field Profile at ', num2str(n), 'th time-step'])
   colorbar
%    rectangle('Position',[-2500e-9,-220e-9,1300e-9,440e-9])
%    rectangle('Position',[-1200e-9,-1200e-9,2400e-9,2400e-9])
%    rectangle('Position',[1200e-9,-720e-9,1300e-9,440e-9])
%    rectangle('Position',[1200e-9,280e-9,1300e-9,440e-9])
   caxis([-2e-10 2e-10])
%    shading interp
   pause(0.01);
%    end
% n
   if(n==1.5e3)
       Ez_2DFigure = Ez;
   end
end
% figure
% pcolor(X*1e9,Y*1e9,Ez_2DFigure);
% % rectangle('Position',[-1500,-90,3000,180])
% rectangle('Position',[-2500e-9,-220e-9,1300e-9,440e-9])
% caxis([-4e-10 4e-10])
% xlabel('X [nm]')
% ylabel('Y [nm]')
% title('Electric Field Profile at The 1500th Time-step')
% shading interp
% colorbar
% rectangle('Position',[-2500,-220,1300,440])
% figure
% set(gcf,'position',[100,100,1500,550])
% plot(log10(abs(Ez_Probe)),'LineWidth',1.5)
% xlabel('Time step','FontSize',12)
% ylabel('log_{10}(|E_{z}|)','FontSize',12)
% title('Logarithm of the Magnitude of the Recorded Electric Field at the Proble Location (x=1000nm,y=0)','FontSize',15)
% grid