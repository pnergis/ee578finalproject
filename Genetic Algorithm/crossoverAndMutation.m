function [populationMatrix,performanceVector] = crossoverAndMutation(populationMatrix,performanceVector,mutationProbability)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
popSize = size(populationMatrix,1);
numofParam = size(populationMatrix,2);

for i=1:2:popSize/2
    crossoverBorder = randi(numofParam-1);
    populationMatrix(popSize/2+i,:) = [populationMatrix(i,1:crossoverBorder),populationMatrix(i+1,crossoverBorder+1:end)];
    for j=1:numofParam
        if(rand<mutationProbability)
            populationMatrix(popSize/2+i,j) = abs(populationMatrix(popSize/2+i,j)-1);
        end
    end
    populationMatrix(popSize/2+i+1,:) = [populationMatrix(i+1,1:crossoverBorder),populationMatrix(i,crossoverBorder+1:end)];
    for j=1:numofParam
        if(rand<mutationProbability)
            populationMatrix(popSize/2+i+1,j) = abs(populationMatrix(popSize/2+i+1,j)-1);
        end
    end
end
performanceVector(popSize/2+1:end) = 0;
end