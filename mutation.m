function U = mutation(p, lu, i, popsize, n, X, F, CR);

% Choose the indices for mutation
indexSet = 1:popsize;
indexSet(i) = [];

% Choose the first Index
temp = floor(rand * (popsize - 1)) + 1;
nouse(1) = indexSet(temp);
indexSet(temp) = [];

% Choose the second index
temp = floor(rand * (popsize - 2)) + 1;
nouse(2) = indexSet(temp);
indexSet(temp) = [];

% Choose the third index
temp = floor(rand * (popsize - 3)) + 1;
nouse(3) = indexSet(temp);

% Mutate
V = p(nouse(1), :) + F.* (p(nouse(2), :) - p(nouse(3), :));

% Handle the elements of the mutant vector which violate the boundary
vioLow = find(V < lu(1, :));
if rand < 0.5
  V(1, vioLow) = lu(1, vioLow);
else
  V(1, vioLow) = 2 .* lu(1, vioLow) - V(1, vioLow);
  vioLowUpper = find(V(1, vioLow) > lu(2, vioLow));
  V(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
end

vioUpper = find(V > lu(2, : ));
if rand < 0.5
  V(1, vioUpper) = lu(2, vioUpper);
else
  V(1, vioUpper) = 2 .* lu(2, vioUpper) - V(1, vioUpper);
  vioUpperLow = find(V(1, vioUpper) < lu(1, vioUpper));
  V(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
end

% Implement the binomial crossover
jRand = floor(rand * n) + 1;
t = rand(1, n) < CR;
t(1, jRand) = 1;
t_ = 1 - t;

U = t .* V + t_ .* X;
