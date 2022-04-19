function efficiency = FDTD_TM(permittivity,enableAnimate,slab,odd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sepehr Eskandari
% EE578: Computational Electromagnetics for Engineers, Spring 2022
% Date: 3/16/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if(~odd)
    [X,Y] = meshgrid(linspace(-xDim/2,xDim/2,Nx),linspace(-yDim/2,yDim/2,Ny));
else

    [X,Y] = meshgrid(linspace(-xDim/2,xDim/2,Nx),linspace(-yDim/2,yDim/2,Ny+1));
end
% Inhomogeneity
er = ones(Ny,Nx);
% Splitter
er(66:185,66:185) = permittivity;
% Upper waveguide
er(151:172,186:end) = 12;
% Bottom waveguide
er(79:100,186:end) = 12;
% Input waveguide
er(115:136,1:65) = 12;

row_observation = 162;
col_observation = Nx-25;
% overwriting the structure with a single slab
if(slab)
    er = ones(Ny,Nx);
    er(115:136,:) = 12;
    row_observation = 126;
    col_observation = Nx-25;
end
% Initialization of fields
Ez = zeros(Ny,Nx);
Dz = zeros(Ny,Nx);
Hx = zeros(Ny+1,Nx);
Bx = zeros(Ny+1,Nx);
Hy = zeros(Ny,Nx+1);
By = zeros(Ny,Nx+1);
% Source parameters
lambdaL = 1534e-9;
lambda0 = 1550e-9;
lambdaU = 1566e-9;
w0 = 2*pi/lambda0;
sigma = 2*lambda0/w0/(lambdaU-lambdaL);
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
Bx_SigmaY = zeros(Ny+1,Nx);
for i=1:L
    Bx_SigmaY(i,:) = maxSigma*(1+(1-i)/L).^gradingOrder;
    Bx_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L).^gradingOrder;
end
Bx1 = (1-Bx_SigmaY*dt/2)./(1+Bx_SigmaY*dt/2);
Bx2 = 1+Bx_SigmaY*dt/2;


% SigmaX evaluated at By location
% By1 and By2 are multiplication factor matrices.
By_SigmaX = zeros(Ny,Nx+1);
for i=1:L
    By_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    By_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
By1 = (1-By_SigmaX*dt/2)./(1+By_SigmaX*dt/2);
By2 = 1+By_SigmaX*dt/2;


% SigmaX evaluated at Hx location
% Hx1 and Hx2 are multiplication factor matrices.
Hx_SigmaX = zeros(Ny+1,Nx);
for i=1:L
    Hx_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Hx_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Hx1 = 1 + Hx_SigmaX*dt/2;
Hx2 = 1 - Hx_SigmaX*dt/2;


% SigmaY evaluated at Hy location
% Hy1 and Hy2 are multiplication factor matrices.
Hy_SigmaY = zeros(Ny,Nx+1);
for i=1:L
    Hy_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Hy_SigmaY(Ny-i+1,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Hy1 = 1 + Hy_SigmaY*dt/2;
Hy2 = 1 - Hy_SigmaY*dt/2;


% SigmaX evaluated at Dz location
% Dz1 and Dz2 are multiplication factor matrices.
Dz_SigmaX = zeros(Ny,Nx);
for i=1:L
    Dz_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Dz_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Dz1 = (1-Dz_SigmaX*dt/2)./(1+Dz_SigmaX*dt/2);
Dz2 = 1+Dz_SigmaX*dt/2;


% SigmaY evaluated at Ez location.
% Ez1 and Ez2 are multiplication factor matrices.
Ez_SigmaY = zeros(Ny,Nx);
for i=1:L
    Ez_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Ez_SigmaY(Ny-i+1,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Ez1 = (1-Ez_SigmaY*dt/2)./(1+Ez_SigmaY*dt/2);
Ez2 = 1+Ez_SigmaY*dt/2;
Ez_observation = zeros(nTotal,1);
Hx_observation = zeros(nTotal,1);
% Update loop
for n=1:nTotal
    Ez_observation(n) = Ez(row_observation,col_observation);
    Hx_observation(n) = Hx(row_observation,col_observation);
    Bx_BeforeUpdate = Bx;
    By_BeforeUpdate = By;
    Dz_BeforeUpdate = Dz;
    Bx(2:Ny,:) = Bx1(2:Ny,:).*Bx(2:Ny,:) - dt./dy./Bx2(2:Ny,:).*(Ez(2:end,:)-Ez(1:(Ny-1),:)); % step 1
    By(:,2:Nx) = By1(:,2:Nx).*By(:,2:Nx) + dt./dx./By2(:,2:Nx).*(Ez(:,2:end)-Ez(:,1:(Nx-1))); % step 2
    Hx(2:Ny,:) = Hx(2:Ny,:) + Hx1(2:Ny,:).*Bx(2:Ny,:) - Hx2(2:Ny,:).*Bx_BeforeUpdate(2:Ny,:); % step 3
    Hy(:,2:Nx) = Hy(:,2:Nx) + Hy1(:,2:Nx).*By(:,2:Nx) - Hy2(:,2:Nx).*By_BeforeUpdate(:,2:Nx); % step 4
    Dz = Dz1.*Dz + dt./Dz2./er.*((Hy(:,2:end)-Hy(:,1:Nx))/dx-(Hx(2:end,:)-Hx(1:Ny,:))/dy); % step 5
    Ez = Ez1.*Ez + 1./Ez2.*(Dz-Dz_BeforeUpdate); % step 6
    if(odd)
        Ez(138:220,L+10) = exp(-0.8636*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma));
        Ez(127:137,L+10) = exp(-0.8636*2.2)/sin(1.0308*2.2)*sin(1.0308e7*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma));
        Ez(126,L+10) = 0;
        Ez(115:125,L+10) = -flipud(Ez(127:137,L+10));
        Ez(32:114,L+10) = -flipud(Ez(138:220,L+10));
    else
        % Hard Source
        Ez(138:220,L+10) = (exp(-1.236*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Ez(127:137,L+10) = (exp(-1.236*2.2)/cos(0.5299*2.2)*cos(5.299e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Ez(126,L+10) = (exp(-1.236*2.2)/cos(0.5299*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Ez(115:125,L+10) = flipud((exp(-1.236*2.2)/cos(0.5299*2.2)*cos(5.299e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
        Ez(32:114,L+10) = flipud((exp(-1.236*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
        % Soft Source
        %         Ez(127:137,L+10) = Ez(127:137,L+10) + (exp(-1.236*2.2)/cos(0.5299*2.2)*cos(5.299e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        %         Ez(126,L+10) = Ez(126,L+10) + (exp(-1.236*2.2)/cos(0.5299*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        %         Ez(115:125,L+10) = Ez(115:125,L+10) + flipud((exp(-1.236*2.2)/cos(0.5299*2.2)*cos(5.299e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
        %         Ez(138:220,L+10) = Ez(138:220,L+10)+ (exp(-1.236*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        %         Ez(32:114,L+10) = Ez(32:114,L+10) + flipud((exp(-1.236*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');

    end

    %%%%%
    if(n>4000 && enableAnimate && mod(n,10)==0 && slab && odd)
        subplot(2,1,1)
        plot(Hx(126,:))
        ylim([-0.3 0.3])
        subplot(2,1,2)
        pcolor(X,Y,Hx);
        title(['Magnetic Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,xDim,440e-9]) % Input Waveguide
        caxis([-0.3 0.3])
        shading interp
        pause(0.01);
    end
    if(n>4000 && enableAnimate && mod(n,10)==0 && ~slab && ~odd)
        pcolor(X,Y,Ez);
        title(['Electric Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,1300e-9,440e-9]) % Input Waveguide
        rectangle('Position',[-xDim/2+1300e-9,-1200e-9,2400e-9,2400e-9]) % Splitter
        rectangle('Position',[-xDim/2+3700e-9,500e-9,xDim-3700e-9,440e-9]) % Output Waveguide - Upper
        rectangle('Position',[-xDim/2+3700e-9,-940e-9,xDim-3700e-9,440e-9])
        caxis([-0.3 0.3])
        shading interp
        pause(0.01);
    end

    if(n>4000 && enableAnimate && mod(n,10)==0 && slab && ~odd)
        subplot(2,1,1)
        plot(Ez(126,:))
        ylim([-0.3 0.3])
        subplot(2,1,2)
        pcolor(X,Y,Ez);
        title(['Electric Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,xDim,440e-9]) % Input Waveguide
        caxis([-0.3 0.3])
        shading interp
        pause(0.01);
    end
    n;
end
if(odd)
    efficiency = max(Hx_observation) / 0.4937
else
    efficiency = max(Ez_observation) / 0.1669; % Use 0.434 for soft sources, and 0.1669 for hard sources
end
end
