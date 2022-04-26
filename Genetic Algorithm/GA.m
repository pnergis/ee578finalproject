% Dont forget to comment s
popSize = 20; % Size of the population, must be a multiply of 4.
numberofIterations = 100;
mutationProbability = 0.1; % Mutation probability.
performanceVector = zeros(popSize,1); % A vector of size popSize.
performanceVectorTM = zeros(popSize,1);
performanceVectorTE = zeros(popSize,1);
% populationMatrix = randi([0 1], popSize,400); % A matrix of size popSize by numberofParams (here 400). Each row corresponds to a specific gene.
populationMatrix = ones(popSize,400); % A matrix of size popSize by numberofParams (here 400). Each row corresponds to a specific gene.
% randomly change 5 bits to 0 of the all ones initial design to minimize randomness
for i=1:popSize
    for j=1:5
        populationMatrix(i,randi([1 400], 1)) = 0;
    end
end
performanceRecordedOverIteration = zeros(numberofIterations,3);
for i=1:numberofIterations
    tic;
    % The first thing to do is to evaluate the performance of our
    % population. However, we only need to evaluate the new genes. So we
    % run the FDTD simulation only for the cases that their performance is
    % zero (has not been calculated yet).
    for j=1:popSize
        if performanceVector(j)==0
            % Each row of the populationMatrix correspond to one sample. To
            % evaluate its performance, first we reshape it to a square of
            % 20 by 20. Then We have to resize this square to a bigger
            % square of size 122 by 122. The reason for that is beacuse our
            % simulation grid is much finer that optimization grid.
            permittivity = 11*imresize(reshape(populationMatrix(j,:),[20,20]),[120,120],'nearest')+1;
            % Now we use the generated profile for the simulation. Our
            % FDTD_TM has the ability to run an animation. Here we passed a
            % 0 as the second argument to disable this ability.
            performanceVectorTE(j) = FDTD_TE(permittivity,0,0,0);
            performanceVectorTM(j) = FDTD_TM(permittivity,0,0,0);
%             if(performanceVectorTE(j) > 1.5*performanceVectorTM(j) || 1.5*performanceVectorTE(j) < performanceVectorTM(j))
%                 performanceVectorTE(j) = 0.7*performanceVectorTE(j);
%                 performanceVectorTM(j) = 0.7*performanceVectorTM(j);
%             end
            performanceVector(j) = (performanceVectorTE(j) + performanceVectorTM(j))/2; 

        end
    end
    % Now we have a generation (populationMatrix) and their performance
    % (performanceVector). No we want to sort our generation in a
    % descending order according to their performance.
    [populationMatrix,performanceVector,performanceVectorTM,performanceVectorTE] = sorterGA(populationMatrix,performanceVector,performanceVectorTM,performanceVectorTE);
    [performanceVector,performanceVectorTM,performanceVectorTE] %#ok<NOPTS> 
    % Now we have a sorted generation. We need to perform crossover and
    % mutation to generate the next generation, and then repeat the for
    % loop.
    performanceRecordedOverIteration(i,1) = performanceVector(1);
    performanceRecordedOverIteration(i,2) = performanceVectorTM(1);
    performanceRecordedOverIteration(i,3) = performanceVectorTE(1);
    [populationMatrix,performanceVector] = crossoverAndMutation(populationMatrix,performanceVector,mutationProbability);
    t2 = toc;
    disp(['Progress: ',num2str(100*i/numberofIterations,'%.2f'),'%'])
    disp(['Time remaining: ', num2str((numberofIterations-i)*t2/60,'%.2f'), ' minutes'])
end
save('GA_DATA_TEeven_TMeven.mat')
