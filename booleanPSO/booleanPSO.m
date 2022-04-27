popSize = 10; 
numberofIterations = 10;
swarm = true([popSize,400]);
for i=1:popSize
    for j=1:5
        swarm(i,randi(400)) = false;
    end
end

velocitySwarm = false([popSize,400]);
performanceVectorCurrent = zeros(popSize,1); % A vector of size popSize used to store the performence of particle in each iteration.
pBest = false([popSize,400]); % A matrix of size popSize by numberofParams (here 400). ith row corresponds to the best position ith particle has been ever explored.
gBest = false([1,400]); % A vector of size 1 by numberofParams. The best overall point ever explored (Best of pBests).
performanceVectorPBest = zeros(popSize,1);
performanceGBest = 0;

% for i=1:numberofIterations
while(performanceGBest<0.85)
    tic;
    % The first thing to do is to evaluate the performance of our swarm
    for j=1:popSize
        % Each row of the populationMatrix correspond to one sample. To
        % evaluate its performance, first we reshape it to a square of
        % 20 by 20. Then We have to resize this square to a bigger
        % square of size 122 by 122. The reason for that is beacuse our
        % simulation grid is much finer that optimization grid.
        permittivity = 11*imresize(reshape(swarm(j,:),[20,20]),[120,120],'nearest')+1;
        % Now we use the generated profile for the simulation. Our
        % FDTD_TM has the ability to run an animation. Here we passed a
        % 0 as the second argument to disable this ability.

        performanceVectorCurrent(j) = FDTD_TM(permittivity,0,0,0);
        if(performanceVectorCurrent(j) > performanceVectorPBest(j))
            performanceVectorPBest(j) = performanceVectorCurrent(j);
            pBest(j,:) = swarm(j,:);
        end
        if(performanceVectorCurrent(j) > performanceGBest)
            performanceGBest = performanceVectorCurrent(j);
            gBest = swarm(j,:);
        end
    end
    % In case that c1 and c2 need to be generated for each particle
    % separately.
    % Now we want to calculate the velocities and update our particles.
    % c1 and c2 are set to be 2.
    for j=1:popSize
        if(rand < 0.5)
            c1 = true;
        else
            c1 = false;
        end
        if(rand < 0.5)
            c2 = true;
        else
            c2 = false;
        end
        if(rand < 0.1)
            w = true;
        else
            w = false;
        end
        velocitySwarm(j,:) = (w&velocitySwarm(j,:)) | (c1&xor(pBest(j,:),swarm(j,:))) | (c1&xor(gBest,swarm(j,:)));
        velocitySwarm(j,:) = speedLimiter(velocitySwarm(j,:),40);
    end
    
    swarm = xor(swarm,velocitySwarm);
    performanceVectorPBest
    performanceGBest
    t2 = toc;
end
save('Data_TM.mat')
