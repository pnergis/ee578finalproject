popSize = 20; % Size of the population, must be a multiply of 4.
numberofIterations = 20;
mutationProbability = 0.1; % Mutation probability.
performanceVector = zeros(popSize,1); % A vector of size popSize.
populationMatrix = randi(2,[popSize,400])-1; % A matrix of size popSize by numberofParams (here 400). Each row corresponds to a specific gene.
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
            permittivity = 11*imresize(reshape(populationMatrix(j,:),[20,20]),[122,122],'nearest')+1;
            % Now we use the generated profile for the simulation. Our
            % FDTD_TM has the ability to run an animation. Here we passed a
            % 0 as the second argument to disable this ability.
            [performanceVector(j),~,~] = FDTD_TM(permittivity,0);
        end
    end
    % Now we have a generation (populationMatrix) and their performance
    % (performanceVector). No we want to sort our generation in a
    % descending order according to their performance.
    [populationMatrix,performanceVector] = sorterGA(populationMatrix,performanceVector);
    performanceVector %#ok<NOPTS> 
    % Now we have a sorted generation. We need to perform crossover and
    % mutation to generate the next generation, and then repeat the for
    % loop.
    [populationMatrix,performanceVector] = crossoverAndMutation(populationMatrix,performanceVector,mutationProbability);
    t2 = toc;
    disp(['Progress: ',num2str(100*i/numberofIterations,'%.2f'),'%'])
    disp(['Time remaining: ', num2str((numberofIterations-i)*t2/60,'%.2f'), ' minutes'])
end
