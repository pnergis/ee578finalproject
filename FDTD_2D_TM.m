
function [Ez, Hy, Hx] = FDTD_2D_TM(pbs)
% Setup

ep0 = 1;                          % permittivity in free space
meu0 = 1;                         % permeability in free space
c0 = 1/sqrt(ep0*meu0);            % speed of light in vacuum
dx = 20e-9;                       % spatial step in x
dy = 20e-9;                       % spatial step in y
Nx = floor(5000e-9/dx);           % number of cells in x
Ny = floor(5000e-9/dy);           % number of cells in y
x = -2500e-9:dx:2500e-9;          % x values in lattice
y = -2500e-9:dy:2500e-9;          % y values in lattice
dt = dx/(sqrt(2)*c0);             % time step
Nt = 1e4;                         % number of steps in time
t = (0:Nt-1)*dt;                  % time
df = 1/(Nt*dt);                   % frequency step
f = 0:df:(Nt-1)*df;               % frequency
lamda = c0./f;                    % wavelength
lamda0 = 1550e-9;                 % center wavelength
BW = 32e-9;                       % bandwidth
lamdaL = lamda0-BW/2;             % lower bound wavelength
lamdaU = lamda0+BW/2;             % upper bound wavelength
w0 = 2*pi*c0/lamda0;              % center frequency
sigma = 2/w0*(lamda0/(lamdaU-lamdaL));

[X,Y] = meshgrid(x(1:end-1),y(1:end-1));

% fields for TM
Ez=zeros(Nx,Ny,Nt);
Dz=zeros(Nx,Ny,Nt);

Hx=zeros(Nx-1,Ny,Nt);
Bx=zeros(Nx-1,Ny,Nt);

By=zeros(Nx,Ny-1,Nt);
Hy=zeros(Nx,Ny-1,Nt);

% sine envelope gaussian
Jz = sin(w0*(t-4*sigma)).*exp(-(t-4*sigma).^2/(sigma^2));

% location of source
src_loc = [125 20];

%% PBS

ep = ones(2*Nx-1, 2*Ny-1)*ep0;
meu = ones(2*Nx-1, 2*Ny-1)*meu0;
wg_size = 440e-9;
wg_dist = 1000e-9; % distance between the center of output wgs
pbs_size = 2400e-9;
epSi = 12;
meuSi = 1;

ep((Ny-wg_size/dy):(Ny+wg_size/dy),1:(Nx-pbs_size/dx)) = epSi;

ep = pixelToEp(pbs, Ny, Nx, dy, dx, ep, epSi, pbs_size);

ep(((Ny+wg_dist/dy)-wg_size/dy):((Ny+wg_dist/dy)+wg_size/dy),(Nx+pbs_size/dx):end) = epSi;
ep(((Ny-wg_dist/dy)-wg_size/dy):((Ny-wg_dist/dy)+wg_size/dy),(Nx+pbs_size/dx):end) = epSi;

meu_Bx = ones(Nx, Nx-1)*meuSi;
meu_By = ones(Ny-1, Ny)*meuSi;
ep_Dz = ones(Nx, Nx);

figure (1)
imagesc(x,y,ep)
xlabel('x')
ylabel('y')
title('\epsilon area of the problem');
set(gca,'YDir','normal')
%% PML

% PML parameters to tune
L = 15; % thickness of the PML
n = 3; % grading of sigma function (given in problem)
R = 1e-9; % power reflection coefficient, typically set to 1e-6

% You must take into account the fact that the nodes Ez, 
% Hx, and Hy are all offset by half a cell in space.

Lx = L*dx;
Ly = L*dy;

x_h1=Lx:-dx/2:0;
x_h2=0:dx/2:Lx;
y_h1=Ly:-dy/2:0;
y_h2=0:dy/2:Ly;

sig_x = zeros(2*Nx, 2*Nx);
sig_y = zeros(2*Ny, 2*Ny);
sig_x_max = -(n+1/4)*(c0/Lx)*log(R);
sig_y_max = -(n+1/4)*(c0/Ly)*log(R);

for i=1:2*Nx-1
    sig_x(i,1:(2*L+1))=abs(sig_x_max*(x_h1(:)/Lx).^n);
    sig_x(i,(2*Nx-1-2*L):(2*Nx-1))=abs(sig_x_max*(x_h2(:)/Lx).^n);
end

for i=1:2*Ny-1
    sig_y(1:(2*L+1),i)=abs(sig_y_max*(y_h1(:)/Ly).^3);
    sig_y((2*Ny-1-2*L):(2*Ny-1),i)=abs(sig_y_max*(y_h2(:)/Ly).^n);
end

sig_y_Ez=zeros(Nx,Ny);
sig_x_Dz=zeros(Nx,Ny);

sig_y_Bx=zeros(Nx-1,Ny);
sig_x_Hx=zeros(Nx-1,Ny);

sig_x_By=zeros(Nx,Ny-1);
sig_y_Hy=zeros(Nx,Ny-1);

for i=1:Nx
    for j=1:Ny-1
        sig_x_By(i,j)=sig_x(2*i-1,2*j);
        sig_y_Hy(i,j)=sig_y(2*i-1,2*j);
        meu_By(i,j)=meu(2*i-1,2*j);
    end
end

for i=1:Nx-1
    for j=1:Ny
       sig_y_Bx(i,j)=sig_y(2*i,2*j-1);
       sig_x_Hx(i,j)=sig_x(2*i,2*j-1);
       meu_Bx(i,j)=meu(2*i,2*j-1);
    end
end


for i=1:Nx
    for j=1:Ny
        sig_x_Dz(i,j)=sig_x(2*i-1,2*j-1);
        ep_Dz(i,j)=ep(2*i-1,2*j-1);
        sig_y_Ez(i,j)=sig_y(2*i-1,2*j-1);
    end
end


%% 2D FDTD

for n=1:Nt-1
    Dz(2:Nx-1,2:Ny-1,n+1)=((1-dt.*sig_x_Dz(2:Nx-1,2:Ny-1)/2)./...
        (1+dt.*sig_x_Dz(2:Nx-1,2:Ny-1)/2)).*Dz(2:Nx-1,2:Ny-1,n)...
          +(dt./((1+dt.*sig_x_Dz(2:Nx-1,2:Ny-1)/2).*...
          ep_Dz(2:Nx-1,2:Ny-1))).*((Hy(2:Nx-1,2:Ny-1,n)-...
          Hy(2:Nx-1,1:Ny-2,n))./dx-(Hx(2:Nx-1,2:Ny-1,n)-...
          Hx(1:Nx-2,2:Ny-1,n))./dy);

    Dz(src_loc(1),src_loc(2),n+1)=Dz(src_loc(1),src_loc(2),n+1)-...
        (dt./ep(src_loc(1),src_loc(2))).*Jz(n);

    Ez(2:Nx-1,2:Ny-1,n+1)=((1-dt.*sig_y_Ez(2:Nx-1,2:Ny-1)/2)./...
        (1+dt.*sig_y_Ez(2:Nx-1,2:Ny-1)/2)).*Ez(2:Nx-1,2:Ny-1,n)+...
        (1./(1+dt.*sig_y_Ez(2:Nx-1,2:Ny-1)/2)).*(Dz(2:Nx-1,2:Ny-1,n+1)-...
        Dz(2:Nx-1,2:Ny-1,n));

    Bx(1:Nx-1,1:Ny,n+1)=((1-dt.*sig_y_Bx(1:Nx-1,1:Ny)/2)./...
        (1+dt.*sig_y_Bx(1:Nx-1,1:Ny)/2)).*Bx(1:Nx-1,1:Ny,n)-...
        (dt./((1+dt.*sig_y_Bx(1:Nx-1,1:Ny)/2).*(meu_Bx(1:Nx-1,1:Ny)))).*...
        ((Ez(2:Nx,1:Ny,n+1)-Ez(1:Nx-1,1:Ny,n+1))./dy); 


    By(1:Nx,1:Ny-1,n+1)=((1-dt.*sig_x_By(1:Nx,1:Ny-1)/2)./...
        (1+dt.*sig_x_By(1:Nx,1:Ny-1)/2)).*By(1:Nx,1:Ny-1,n)+...
        (dt./((1+dt.*sig_x_By(1:Nx,1:Ny-1)/2).*...
        (meu_By(1:Nx,1:Ny-1)))).*...
        ((Ez(1:Nx,2:Ny,n+1)-Ez(1:Nx,1:Ny-1,n+1))./dx);    


    Hx(1:Nx-1,1:Ny,n+1)=Hx(1:Nx-1,1:Ny,n)+...
        (1+dt.*sig_x_Hx(1:Nx-1,1:Ny)/2).*Bx(1:Nx-1,1:Ny,n+1)-...
        (1-dt.*sig_x_Hx(1:Nx-1,1:Ny)/2).*Bx(1:Nx-1,1:Ny,n);

    Hy(1:Nx,1:Ny-1,n+1)= Hy(1:Nx,1:Ny-1,n)+...
        (1+dt.*sig_y_Hy(1:Nx,1:Ny-1)/2).*By(1:Nx,1:Ny-1,n+1)-...
        (1-dt.*sig_y_Hy(1:Nx,1:Ny-1)/2).*By(1:Nx,1:Ny-1,n);
    
    if(n>3000)
       pcolor(X,Y,Ez(:,:,n));
       title(['Electric Field Profile at ', num2str(n), 'th time-step'])
       colorbar
       rectangle('Position',[-2500e-9,-220e-9,1300e-9,440e-9])
       rectangle('Position',[-1200e-9,-1200e-9,2400e-9,2400e-9])
       rectangle('Position',[1200e-9,280e-9,1300e-9,440e-9])
       rectangle('Position',[1200e-9,-720e-9,1300e-9,440e-9])
       caxis([-4e-10 4e-10])
       shading interp
       pause(0.02);
   end
end
end