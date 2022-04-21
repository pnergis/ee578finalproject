function [populationMatrix,performanceVector,performanceVectorTM,performanceVectorTE] = sorterGA(populationMatrix,performanceVector,performanceVectorTM,performanceVectorTE);
%Summary of this function goes here
%   populationMatrix is a popSize by numofParam matrix, where popSize is
%   the population size in genetic algorithm, and numofParam is the number
%   of parameters to be optimized. 
%   performanceVector is a vector of size popSize.
%   Our goal here is to sort the rows of populationMatrix in a descending
%   order according to their correspoding performance in performanceVector.
[performanceVector,I] = sort(performanceVector,'descend');
performanceVectorTM = performanceVectorTM(I);
performanceVectorTE = performanceVectorTE(I);
populationMatrix = populationMatrix(I,:);
end
