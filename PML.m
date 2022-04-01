function []

% PML parameters to tune
L = 10; % thickness of the PML
n = 3; % grading of sigma function (given in problem)
R = 1e-6; % power reflection coefficient, typically set to 1e-6

% You must take into account the fact that the nodes Ez, 
% Hx, and Hy are all offset by half a cell in space.

Lx = L*dx;
Ly = L*dy;

x_h1=Lx:-dx/2:0;
x_h2=0:dx/2:Lx;
y_h1=Ly:-dy/2:0;
y_h2=0:dy/2:Ly;

ep = ones(2*Ny-1, 2*Nx-1)*ep0;
meu = ones(2*Ny-1, 2*Nx-1)*meu0;
meu_Bx = ones(Nx, Nx-1);
meu_By = ones(Ny-1, Ny);
ep_Dz = ones(Nx, Nx);

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
sig_x_Hxtilda=zeros(Nx-1,Ny);

sig_x_By=zeros(Nx,Ny-1);
sig_y_Hytilda=zeros(Nx,Ny-1);

for i=1:Nx
    for j=1:Nx-1
        sig_x_By(i,j)=sig_x(2*i-1,2*j);
        sig_y_Hytilda(i,j)=sig_y(2*i-1,2*j);
        meu_By(i,j)=meu(2*i-1,2*j);
    end
end

for i=1:Nx-1
    for j=1:Nx
       sig_y_Bx(i,j)=sig_y(2*i,2*j-1);
       sig_x_Hxtilda(i,j)=sig_x(2*i,2*j-1);
       meu_Bx(i,j)=meu(2*i,2*j-1);
    end
end


for i=1:Nx
    for j=1:Nx
        sig_x_Dz(i,j)=sig_x(2*i-1,2*j-1);
        ep_Dz(i,j)=ep(2*i-1,2*j-1);
        sig_y_Ez(i,j)=sig_y(2*i-1,2*j-1);
    end
end
