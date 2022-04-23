% Typical use: FDTD_TE(permittivity,0,0,0) -> Use even mode. Odd mode need
% to be double checked.
% For soft source implementation uncomment lines 202 to 206 and modify the
% line 262 (Refer to its comment).
function efficiency = FDTD_TE(permittivity,enableAnimate,slab,odd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sepehr Eskandari
% EE578: Computational Electromagnetics for Engineers, Spring 2022
% Date: 3/16/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to be change in case of increasing the size of domain:
% row_observation, col_observation
% Domain Size
xDim = 7000e-9;
yDim = 5000e-9;
% Grid step size
dx = 20e-9;
dy = 20e-9;
% Time stepd
dt = dx/sqrt(2);
% Total number of time steps
nTotal = 1e4;
% Number of grids in each direction
Nx = floor(xDim/dx);
Ny = floor(yDim/dy);
if(~odd)
    [X,Y] = meshgrid(linspace(-xDim/2,xDim/2,Nx+1),linspace(-yDim/2+dy/2,yDim/2-dy/2,Ny+1));
elseif(odd)
    [X,Y] = meshgrid(linspace(-xDim/2,xDim/2,Nx+1),linspace(-yDim/2+dy/2,yDim/2-dy/2,Ny));
end
% Inhomogeneity

er_Dx = ones(Ny,Nx+1);
er_Dy = ones(Ny+1,Nx);

% er_Dx
% We've divided our splitter to 120x120 cells, and we have permittivity at
% the center of these cells. In order to calculate permittivity at the
% cells edge we have to take the average of neighbors cells.

% Input waveguide
er_Dx(115:136,1:65) = 12;
% Upper output waveguide
er_Dx(151:172,187:end) = 12;
% Bottom output waveguide
er_Dx(79:100,187:end) = 12;
% Splitter
temp1 = [er_Dx(66:185,65),permittivity];
temp2 = [permittivity,er_Dx(66:185,187)];
er_Dx(66:185,66:186) = (temp1+temp2)/2;


% er_Dy
% Input waveguide
er_Dy(116:136,1:65) = 12;
er_Dy(115,1:65) = 6.5;
er_Dy(137,1:65) = 6.5;
% Upper output waveguide
er_Dy(152:172,186:end) = 12;
er_Dy(151,186:end) = 6.5;
er_Dy(173,186:end) = 6.5;
% Bottom output waveguide
er_Dy(80:100,186:end) = 12;
er_Dy(79,186:end) = 6.5;
er_Dy(101,186:end) = 6.5;
% Splitter
temp1 = [permittivity; ones(1,120)];
temp2 = [ones(1,120); permittivity];
er_Dy(66:186,66:185) = (temp1+temp2)/2; % Splitter

row_observation = 90;
col_observation = Nx-25;
% overwriting the structure with a single slab
if(slab)
    er_Dx = ones(Ny,Nx+1);
    er_Dy = ones(Ny+1,Nx);
    er_Dx(115:136,:) = 12;
    er_Dy(116:136,:) = 12;
    er_Dy(115,:) = 6.5;
    er_Dy(137,:) = 6.5;
    row_observation = 126;
    col_observation = Nx-25;
end

% Initialization of fields
Hz = zeros(Ny+1,Nx+1);
Bz = zeros(Ny+1,Nx+1);
Ex = zeros(Ny,Nx+1);
Dx = zeros(Ny,Nx+1);
Ey = zeros(Ny+1,Nx);
Dy = zeros(Ny+1,Nx);
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

% SigmaY evaluated at Dx location
% Dx1 and Dx2 are multiplication factor matrices.
Dx_SigmaY = zeros(Ny,Nx+1);
for i=1:L
    Dx_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Dx_SigmaY(Ny-i+1,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Dx1 = (1-Dx_SigmaY*dt/2)./(1+Dx_SigmaY*dt/2);
Dx2 = 1+Dx_SigmaY*dt/2;


% SigmaX evaluated at Dy location
% Dy1 and Dy2 are multiplication factor matrices.
Dy_SigmaX = zeros(Ny+1,Nx);
for i=1:L
    Dy_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    Dy_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
end
Dy1 = (1-Dy_SigmaX*dt/2)./(1+Dy_SigmaX*dt/2);
Dy2 = 1+Dy_SigmaX*dt/2;


% SigmaX evaluated at Ex location
% Ex1 and Ex2 are multiplication factor matrices.
Ex_SigmaX = zeros(Ny,Nx+1);
for i=1:L
    Ex_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Ex_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Ex1 = 1 + Ex_SigmaX*dt/2;
Ex2 = 1 - Ex_SigmaX*dt/2;


% SigmaY evaluated at Ey location
% Ey1 and Ey2 are multiplication factor matrices.
Ey_SigmaY = zeros(Ny+1,Nx);
for i=1:L
    Ey_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Ey_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Ey1 = 1 + Ey_SigmaY*dt/2;
Ey2 = 1 - Ey_SigmaY*dt/2;


% SigmaX evaluated at Bz location
% Bz1 and Bz2 are multiplication factor matrices.
Bz_SigmaX = zeros(Ny+1,Nx+1);
for i=1:L
    Bz_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Bz_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Bz1 = (1-Bz_SigmaX*dt/2)./(1+Bz_SigmaX*dt/2);
Bz2 = 1+Bz_SigmaX*dt/2;


% SigmaY evaluated at Hz location.
% Hz1 and Hz2 are multiplication factor matrices.
Hz_SigmaY = zeros(Ny+1,Nx+1);
for i=1:L
    Hz_SigmaY(i,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
    Hz_SigmaY(Ny+1-i+1,:) = maxSigma*(1+(1-i)/L)^gradingOrder;
end
Hz1 = (1-Hz_SigmaY*dt/2)./(1+Hz_SigmaY*dt/2);
Hz2 = 1+Hz_SigmaY*dt/2;

Ex_observation = zeros(nTotal,1);
Hz_observation = zeros(nTotal,1);
% Update loop
for n=1:nTotal
    Ex_observation(n) = Ex(row_observation,col_observation);
    Hz_observation(n) = Hz(row_observation,col_observation);
    Dx_BeforeUpdate = Dx;
    Dy_BeforeUpdate = Dy;
    Bz_BeforeUpdate = Bz;
    Dx(:,2:Nx) = Dx1(:,2:Nx).*Dx(:,2:Nx) + dt./dy./Dx2(:,2:Nx)./er_Dx(:,2:Nx).*(Hz(2:end,2:Nx)-Hz(1:Ny,2:Nx)); % step 1
    Dy(2:Ny,:) = Dy1(2:Ny,:).*Dy(2:Ny,:) - dt./dx./Dy2(2:Ny,:)./er_Dy(2:Ny,:).*(Hz(2:Ny,2:end)-Hz(2:Ny,1:Nx)); % step 2
    Ex(:,2:Nx) = Ex(:,2:Nx) + Ex1(:,2:Nx).*Dx(:,2:Nx) - Ex2(:,2:Nx).*Dx_BeforeUpdate(:,2:Nx); % step 3
    Ey(2:Ny,:) = Ey(2:Ny,:) + Ey1(2:Ny,:).*Dy(2:Ny,:) - Ey2(2:Ny,:).*Dy_BeforeUpdate(2:Ny,:); % step 4
    Bz(2:Ny,2:Nx) = Bz1(2:Ny,2:Nx).*Bz(2:Ny,2:Nx) + dt./Bz2(2:Ny,2:Nx).*((Ex(2:end,2:Nx)-Ex(1:(Ny-1),2:Nx))./dy-(Ey(2:Ny,2:Nx)-Ey(2:Ny,1:(Nx-1)))./dx); % step 5
    Hz(2:Ny,2:Nx) = Hz1(2:Ny,2:Nx).*Hz(2:Ny,2:Nx) + 1./Hz2(2:Ny,2:Nx).*(Bz(2:Ny,2:Nx)-Bz_BeforeUpdate(2:Ny,2:Nx)); % step 6
    if(odd)
        % Hard Sources
        Hz(138:220,L+10) = (exp(-3.5964*1e6*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Hz(127:137,L+10) = (exp(-0.35964*2.2)/sin(1.2955*2.2)*sin(1.2955e7*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Hz(126,L+10) = 0;
        Hz(115:125,L+10) = -flipud(Hz(127:137,L+10));
        Hz(32:114,L+10) = -flipud(Hz(138:220,L+10));
    else
        %         Hard Sources
        Hz(138:220,L+10) = (exp(-1.153*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Hz(127:137,L+10) = (exp(-1.153*2.2)/cos(0.6913*2.2)*cos(6.913e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Hz(126,L+10) = (exp(-1.153*2.2)/cos(0.6913*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        Hz(115:125,L+10) = flipud((exp(-1.153*2.2)/cos(0.6913*2.2)*cos(6.913e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
        Hz(32:114,L+10) = flipud((exp(-1.153*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
        % Soft Sources
        %         Hz(138:220,L+10) = Hz(138:220,L+10) + (exp(-1.153*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        %         Hz(127:137,L+10) = Hz(127:137,L+10) + (exp(-1.153*2.2)/cos(0.6913*2.2)*cos(6.913e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        %         Hz(126,L+10) = Hz(126,L+10) + (exp(-1.153*2.2)/cos(0.6913*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
        %         Hz(115:125,L+10) = Hz(115:125,L+10) + flipud((exp(-1.153*2.2)/cos(0.6913*2.2)*cos(6.913e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
        %         Hz(32:114,L+10) = Hz(32:114,L+10) + flipud((exp(-1.153*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
    end
    if(n>4000 && mod(n,10)==0 && enableAnimate && ~slab && odd)
        pcolor(X,Y,Ex);
        title(['Magnetic Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,1300e-9,440e-9]) % Input Waveguide
        rectangle('Position',[-xDim/2+1300e-9,-1200e-9,2400e-9,2400e-9]) % Splitter
        rectangle('Position',[-xDim/2+3700e-9,500e-9,xDim-3700e-9,440e-9]) % Output Waveguide - Upper
        rectangle('Position',[-xDim/2+3700e-9,-940e-9,xDim-3700e-9,440e-9])
        caxis([-1 1])
        shading interp
        pause(0.01);
    end
    if(n>4000 && mod(n,10)==0 && enableAnimate && slab && odd)
        subplot(2,1,1)
        plot(Ex(126,:))
        ylim([-1.7 1.7])
        subplot(2,1,2)
        pcolor(X,Y,Ex);
        title(['Electric Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,xDim,440e-9]) % Input Waveguide
        caxis([-1 1])
        shading interp
        pause(0.01);
    end
    if(n>4000 && mod(n,10)==0 && enableAnimate && ~slab && ~odd)
        pcolor(X,Y,Hz);
        title(['Magnetic Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,1300e-9,440e-9]) % Input Waveguide
        rectangle('Position',[-xDim/2+1300e-9,-1200e-9,2400e-9,2400e-9]) % Splitter
        rectangle('Position',[-xDim/2+3700e-9,500e-9,xDim-3700e-9,440e-9]) % Output Waveguide - Upper
        rectangle('Position',[-xDim/2+3700e-9,-940e-9,xDim-3700e-9,440e-9])
        caxis([-1 1])
        shading interp
        pause(0.01);
    end
    if(n>4000 && mod(n,10)==0 && enableAnimate && slab && ~odd)
        subplot(2,1,1)
        plot(Hz(126,:))
        ylim([-1.7 1.7])
        subplot(2,1,2)
        pcolor(X,Y,Hz);
        title(['Magnetic Field Profile at ', num2str(n), 'th time-step'])
        colorbar
        rectangle('Position',[-xDim/2,-220e-9,xDim,440e-9]) % Input Waveguide
        caxis([-1 1])
        shading interp
        pause(0.01);
    end
end
if(odd)
    efficiency = max(Ex_observation)/0.4156;
else
    efficiency = max(Hz_observation)/1.5809; % Use 4.4143 for soft sources, and 1.5809 for hard source.
end
end
