function [efficiency, er, Hz, Ey, Ex] = FDTD_2D_TE(pbs)
    % Domain Size
    xDim = 5000e-9;
    yDim = 5000e-9;
    % Grid step size
    dx = 20e-9;
    dy = 20e-9;
    % Time step
    dt = dx/sqrt(2);
    % Total number of time steps
    nTotal = 1e4;
    % Number of grids in each direction
    Nx = floor(xDim/dx);
    Ny = floor(yDim/dy);    
    % PBS pixel to epsilon grid
    er_Dx = ones(Ny+1,Nx);
    er_Dy = ones(Ny,Nx+1);
    er = ones(Ny+1,Nx+1);
    erSi = 12;
    wg_size = 440e-9;
    wg_dist = 1000e-9; % distance between the center of output wgs
    pbs_size = 2400e-9;
    %   er(151:173,186:end) = 12;
    % 	er(79:101,186:end) = 12;
    er((Ny/2-wg_size/dy/2)+1:(Ny/2+wg_size/dy/2),...
        1:(Nx/2-pbs_size/dx/2)) = erSi;
    er(((Ny/2+wg_dist/dy/2)-wg_size/dy/2):((Ny/2+wg_dist/dy/2)+wg_size/dy/2),...
        (Nx/2+pbs_size/dx/2)+1:end) = erSi;
    er(((Ny/2-wg_dist/dy/2)-wg_size/dy/2):((Ny/2-wg_dist/dy/2)+wg_size/dy/2),...
        (Nx/2+pbs_size/dx/2)+1:end) = erSi;
    er = pixelToEr(pbs, er, erSi, Ny, Nx, dy, dx, pbs_size); % pbs
    
    er_Dy = er(1:end-1,:);
    er_Dx = er(:,1:end-1);
    % Initialization of fields
    Hz = zeros(Ny,Nx);
    Bz = zeros(Ny,Nx);
    Ex = zeros(Ny+1,Nx);
    Dx = zeros(Ny+1,Nx);
    Ey = zeros(Ny,Nx+1);
    Dy = zeros(Ny,Nx+1);
    % Source parameters
    lambdaL = 1534e-9;
    lambda0 = 1550e-9;
    lambdaU = 1566e-9;
    w0 = 2*pi/lambda0;
    sigma = 2*lambda0/w0/(lambdaU-lambdaL);
    % Source and probe locations
    xSource = -1850e-9;
    ySource = 0;
    x_observation = 2020e-9;
    y_observation = -500e-9;
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
    % Dy1 and Dy2 are multiplication factor matrices.
    Dy_SigmaX = zeros(Ny,Nx+1);
    for i=1:L
        Dy_SigmaX(:,i) = maxSigma*(1+(1-i)/L)^gradingOrder;
        Dy_SigmaX(:,Nx+1-i+1) = maxSigma*(1+(1-i)/L)^gradingOrder;
    end
    Dy1 = (1-Dy_SigmaX*dt/2)./(1+Dy_SigmaX*dt/2);
    Dy2 = 1+Dy_SigmaX*dt/2;


    % SigmaX evaluated at Ex location
    % Ex1 and Ex2 are multiplication factor matrices.
    Ex_SigmaX = zeros(Ny+1,Nx);
    for i=1:L
        Ex_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
        Ex_SigmaX(:,Nx+1-i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    end
    Ex1 = 1 + Ex_SigmaX*dt/2;
    Ex2 = 1 - Ex_SigmaX*dt/2;


    % SigmaY evaluated at Ey location
    % Ey1 and Ey2 are multiplication factor matrices.
    Ey_SigmaY = zeros(Ny,Nx+1);
    for i=1:L
        Ey_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
        Ey_SigmaY(Ny+1-i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    end
    Ey1 = 1 + Ey_SigmaY*dt/2;
    Ey2 = 1 - Ey_SigmaY*dt/2;


    % SigmaX evaluated at Bz location
    % Bz1 and Bz2 are multiplication factor matrices.
    Bz_SigmaX = zeros(Ny,Nx);
    for i=1:L
        Bz_SigmaX(:,i) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
        Bz_SigmaX(:,Nx-i+1) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    end
    Bz1 = (1-Bz_SigmaX*dt/2)./(1+Bz_SigmaX*dt/2);
    Bz2 = 1+Bz_SigmaX*dt/2;


    % SigmaY evaluated at Hz location.
    % Hz1 and Hz2 are multiplication factor matrices.
    Hz_SigmaY = zeros(Ny,Nx);
    for i=1:L
       Hz_SigmaY(i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
       Hz_SigmaY(Ny+1-i,:) = maxSigma*(1+(1/2-i)/L)^gradingOrder;
    end
    Hz1 = (1-Hz_SigmaY*dt/2)./(1+Hz_SigmaY*dt/2);
    Hz2 = 1+Hz_SigmaY*dt/2;


    Hz_Probe = zeros(nTotal,1);
    Hz_Excitation = zeros(nTotal,1);
    % Update loop
    J = zeros(nTotal,1);
    for n=1:nTotal
       Dx_BeforeUpdate = Dx;
       Dy_BeforeUpdate = Dy;
       Bz_BeforeUpdate = Bz;
       Dx(2:Ny,:) = Dx1(2:Ny,:).*Dx(2:Ny,:) + dt./dy./Dx2(2:Ny,:)./er_Dx(2:Ny,:).*(Hz(2:Ny,:)-Hz(1:(Ny-1),:)); 
       Dy(:,2:Nx) = Dy1(:,2:Nx).*Dy(:,2:Nx) - dt./dx./Dy2(:,2:Nx)./er_Dy(:,2:Nx).*(Hz(:,2:Nx)-Hz(:,1:(Nx-1))); 
       Ex = Ex + Ex1.*Dx - Ex2.*Dx_BeforeUpdate; 
       Ey = Ey + Ey1.*Dy - Ey2.*Dy_BeforeUpdate; 
       Bz = Bz1.*Bz + dt./Bz2.*((Ex(2:end,:)-Ex(1:Ny,:))./dy-(Ey(:,2:end)-Ey(:,1:Nx))./dx); 
       Hz = Hz1.*Hz + 1./Hz2.*(Bz-Bz_BeforeUpdate);
       %%%%% Fundamental Even TE Mode
       Hz(127:137,L+10) = Hz(127:137,L+10) + (exp(-1.153*2.2)/cos(0.6913*2.2)*cos(6.913e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
       Hz(126,L+10) = Hz(126,L+10) + (exp(-1.153*2.2)/cos(0.6913*2.2)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
       Hz(115:125,L+10) = Hz(115:125,L+10) + flipud((exp(-1.153*2.2)/cos(0.6913*2.2)*cos(6.913e6*(1:11)*dy)*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))');
       Hz(138:220,L+10) = Hz(138:220,L+10) + (exp(-1.153*1e7*dy*(12:94))*exp(-(n*dt-4*sigma)^2/(sigma^2))*sin(w0*(n*dt-4*sigma)))';
       Hz(32:114,L+10) = Hz(32:114,L+10) + flipud(Hz(138:220,L+10));
       
       %%%%%
       Hz_Probe(n) = Hz(row_observation,col_observation);
       Hz_Excitation(n) = Hz(126,L+10);
%        if(n>3000)
%            pcolor(Hz);
%            title(['Magnetic Field Profile at ', num2str(n), 'th time-step'])
%            colorbar
%            rectangle('Position',[-2500e-9,-220e-9,1300e-9,440e-9])
%            rectangle('Position',[-1200e-9,-1200e-9,2400e-9,2400e-9])
%            rectangle('Position',[1200e-9,-720e-9,1300e-9,440e-9])
%            rectangle('Position',[1200e-9,280e-9,1300e-9,440e-9])
%            caxis([-0.01 0.01])
%            shading interp
%            pause(0.01);
%        end
    end
    efficiency = 100*(max(Hz_Probe) / max(Hz_Excitation));
end
