% Binary Particle Swarm Optimization
% only input how many particles

function [FOM_TE, FOM_TM, pbs] = binPSO(N)
    
    % i refers to particle index and j refers pixel in the
    % reshaped pbs 20x20 design

    % initialize the swarm Xi, the positions of the particles are randomly
    % initialized within the hypercube. Elements of Xi are randomly
    % selected from binary values 0 and 1
    X = randi([0 1], N, 400);
    
    % 2 element row for TE FOM and TM FOM values then 400 element pbs]
    pbest = zeros(N,402); 
    pGbest = zeros(1,402);
    
    d1 = zeros(N,400);
    d0 = zeros(N,400);
    d1G = zeros(N,400);
    d0G = zeros(N,400);
    
    v1 = zeros(N,400);
    v0 = zeros(N,400);
    vC = zeros(N,400); 
    r = zeros(N,400);
    
    % initialize constants
    c1 = 2; c2 = 2;
    w=0.5;
    

    % preallocate for 1000 iterations, this will save the TE and TM
    % efficiency for each particles FDTDs at each iteration of the while
    % loop
    FOM_TE = zeros(N,1000);
    FOM_TM = zeros(N,1000);
    
    % Begin the conditional while loop
    t = 1;
    % begin optimization loop which ends once TE and TM have minimum 85%
    % field transmission efficiency each in their respective output waveguides
    while (pGbest(1)<0.70 || pGbest(2)<0.70)
        r1 = randi([0 1]); r2 = randi([0 1]); % updated each iteration
        
        % Evaluate the performace FOMs of each particle using its current
        % positions Xi(t)
        
        % run FDTD for new FOMs
        % keep track of each particles improvement over iterations
        % to be plotted for report
        for i=1:N
            pbs = 11*imresize(reshape(X(i,:),[20,20]),[120,120],'nearest')+1;
            FOM_TM(i,t) = FDTD_2D_TM(pbs,0,0,0);
            FOM_TE(i,t) = FDTD_2D_TE(pbs,0,0,0);
        end
        
        % Compare the performance of each individual to its best performance
        % so far
        for i=1:N
            if FOM_TE(i,t) > pbest(i,1) && FOM_TM(i,t) > pbest(i,2)
                pbest(i,1) = FOM_TE(i,t);
                pbest(i,2) = FOM_TM(i,t);
                pbest(i,3:end) = X(i,:);
            end
        end
        
        % Compare the performance of each particle  to the global best
        % particle
        for i=1:N
            if FOM_TE(i,t) > pGbest(1) && FOM_TM(i,t) > pGbest(2)
                pGbest(1) = FOM_TE(i,t);
                pGbest(2) = FOM_TM(i,t);
                pGbest(3:end) = X(i,:);
            end
        end
        
        % Change the velocity of the particle, v0(i), v1(i) according to eq
        % (6,7) from the paper
        
        for i =1:N
            for j = 1:400
                % update individual particle best update eq constants 
                if pbest(i,j+2) == 1
                    d1(i,j) = c1*r1;
                    d0(i,j) = -c1*r1;
                else
                    d1(i,j) = -c1*r1;
                    d0(i,j) = c1*r1;
                end

                % update global particle best update eq constants
                if pGbest(j+2) == 1
                    d1G(i,j) = c2*r2;
                    d0G(i,j) = -c2*r2;
                else
                    d1G(i,j) = -c1*r1;
                    d0G(i,j) = c1*r1;
                end

                v1(i,j) = w*v1(i,j) + d1(i,j) + d1G(i,j);
                v0(i,j) = w*v0(i,j) + d0(i,j) + d0G(i,j);

                % Calculate the velocity of change of the bits, Vc(i) as in (5)
                % set up bit flip
                if X(i,j) == 0
                    vC(i,j) = v1(i,j);
                else
                    vC(i,j) = v0(i,j);
                end
                
                % Normalize vC
                vC(i,j) = 1/(1+exp(-vC(i,j)));
                
                % Generate the random variable rij in the range (0,1). Move
                % each particle to a new position using eq (8)
                % to flip or not to flip to the 2's complement
                r(i,j) = randi([0 1]);
                if r(i,j) < vC(i,j)
                    if X(i,j)==0
                        X(i,j) = 1;
                    else
                        X(i,j) = 0;
                    end
                end
            end
            % go to step 2, and repeat until converagence. Aka continue the
            % while loop
        end
        disp([t pGbest(1) pGbest(2)])
        figure (1)
        pcolor(pbs)
        xlabel('x')
        ylabel('y') 
        title('\epsilon');
        set(gca,'YDir','normal')
        t=t+1;
    end
end