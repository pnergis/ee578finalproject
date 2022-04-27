function velocityVector = speedLimiter(velocityVector,Vmax)
V = sum(velocityVector);
if(V>Vmax)
    nonZerosElementsIndices = find(velocityVector); % Obtaining the nonzero elements in velocity vector.
    p = randperm(length(nonZerosElementsIndeces)); % Creating a randomely permutated vector so we can make zero additional elements randomely.
    for i=1:(V-Vmax)
        velocityVector(nonZerosElementsIndices(p(i))) = 0;
    end
end
end
