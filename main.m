%**************************************************************************************************
%Author: Yong Wang
%Last Edited: 01/09/2008
%Email: ywang@csu.edu.cn; wangyong1226@gmail.com
%Reference:          A Hybrid Multi-Swarm Particle Swarm Optimization to Solve
%                                             Constrained Optimization Problems
%                              Frontiers of Computer Science in China, 2009, 3(1):38-52.
%**************************************************************************************************

% profile on;

clc;
clear all;
tic;

format long;
format compact;

% choose which problem to be tested
problemSet = [1 : 19 21 23 24];


    problem = problemSet(problemIndex)

    % Record the best results
    bestResults = [];
    
    % Record the feasibility proportion of the final swarm
    feasiPro = [];

    time = 1;
    
    % The total number of runs
    totalTime = 1;

    while time <= totalTime

        switch problem

            case 1

                % min_var and max_var define the lower and upper bounds of the particle, respectively
                min_var = [0 0 0 0 0 0 0 0 0 0 0 0 0]; max_var = [1 1 1 1 1 1 1 1 1 100 100 100 1];
                % n is the dimension of the particle
                n = 13;
                A = []; %\\\\ - 15\

        end

        rand('seed', sum(100 * clock));

        % The size of the main swarm
        popsize = 60;
        % The size of each sub-swarm
        subpopsize = 8;

        % The tolerance value for equality constraints
        delta = 0.0001;

        totalGen = 300000/ 2 / popsize;

        % Initialize the main swarm
        p = ones(popsize, 1) * min_var + rand(popsize, n) .* (ones(popsize, 1) * (max_var - min_var));

        % Evaluate the main swatm
        fit = fitness(p, problem, delta, A);

        % Determine the pbest
        pbest = p;
        fitPbest = fit;

        gen = 1;

        while gen <= totalGen

            %................ global search model ..............%

            for i = 1:popsize

                % Use DE to update the pbest at each generation
                % In this case, each pbest is considered as a target vector
                % X of DE, U denotes the trial vector
                X = pbest(i, :);
                fitX = fitPbest(i, :);

                % The parameter settings for DE
                F = 0.7;
                CR = 1.0;

                % Implement DE to generate the trial vector
                U = mutation(pbest, [min_var; max_var], i, popsize, n, X, F, CR);

                % Evaluate the trial vector
                fitU = fitness(U, problem, delta, A);

                % Choose the better one (the feasibility criterion) between
                % the trial vector and the target vector to enter the next population
                if (fitX(1, 2) == 0 & fitU(1, 2) == 0)

                    if fitX(1, 1) > fitU(1, 1)
                        pbest(i, :) = U;
                        fitPbest(i, :) = fitU;
                    end

                else

                    if fitX(1, 2) > fitU(1, 2)
                        pbest(i, :) = U;
                        fitPbest(i, :) = fitU;
                    end

                end

            end

            %................ local search model ..............%

            % 1. Sort the swarm in ascending order by the degree of constraint
            % violation
            [sortVal, sortIndex] = sort(fit(:, 2));
            p = p(sortIndex, :);
            fit = fit(sortIndex, :);

            % 2. Sort the feasible particles in ascending order by their objective function
            % values
            findIndex = find(fit(:, 2) == 0);
            if length(findIndex)>0
                [sortVal, sortIndex] = sort(fit(findIndex, 1));
                p(1 : length(findIndex), :) = p(sortIndex, :);
                fit(1 : length(findIndex), :) = fit(sortIndex, :);
            end

            p_ = [];
            pbest_ = [];
            fitPbest_ = [];

            k = 0;
            while k<  floor(popsize / subpopsize)

                % The first particle is the local best of the sub-swarm to which it belong
                lbest = pbest(1, :);

                psize = size(p, 1);

                % Compute the distance from the particles in the population
                % to the local best
                [neighborVal, neighborIndex] = sort(sum((ones(psize, 1) * p(1, :) - p).^2, 2));

                % Several particles with the largest Euclidean distance from
                % the local best are assigned as the other members of each
                % sub-swarm
                subpop = [p(1, :); p(neighborIndex(psize - subpopsize + 2 : psize), :)];
                fitSubpop = [fit(1, :); fit(neighborIndex(psize - subpopsize + 2 : psize), :)];

                % Obtain the pbests of the particles in the sub-swarm
                pbest_s = [pbest(1, :); pbest(neighborIndex(psize - subpopsize + 2 : psize), :)];
                fitPbest_s = [fitPbest(1, :); fitPbest(neighborIndex(psize - subpopsize + 2 : psize), :)];

                p(neighborIndex(psize - subpopsize + 2 : psize), :) = [];
                p(1, :) = [];
                fit(neighborIndex(psize - subpopsize + 2 : psize), :) = [];
                fit(1, :) = [];

                pbest(neighborIndex(psize - subpopsize + 2 : psize), :) = [];
                pbest(1, :) = [];
                fitPbest(neighborIndex(psize - subpopsize + 2 : psize), :) = [];
                fitPbest(1, :) = [];

                % The local best of each sub-swarm undergoes
                % the velocity and position updates with a probability of
                % 0.85
                if rand < 0.85

                    for i = 1 : subpopsize

                        % Initialize the velocity for each particle in
                        % the sub-swarm
                        vel = zeros(subpopsize, n);

                        % Update the velocity of each particle
                        vel(i, :) = abs(randn(1, n)) .* (pbest_s(i, :) - subpop(i, :)) + abs(randn(1, n)) .* (lbest - subpop(i, :));

                        temp_s = subpop(i, :);

                        % Update the position of each particle
                        subpop(i, :) = subpop(i, :) + vel(i, :);

                        % Handle the elements of the particle which violate the boundary
                        tag1 = find(subpop(i, :) < min_var);
                        tag2 = find(subpop(i, :) > max_var);

                        if length(tag1)>0
                            subpop(i, tag1) = (temp_s(1, tag1) + min_var(1, tag1)) .* 0.5;
                        end

                        if length(tag2)>0
                            subpop(i, tag2) = (temp_s(1, tag2) + max_var(1, tag2)) .* 0.5;
                        end

                    end

                else

                    for i = 2 : subpopsize

                        % Initialize the velocity for each particle in
                        % the sub-swarm
                        vel = zeros(subpopsize - 1, n);

                        % Update the velocity of each particle
                        vel(i, :) = abs(randn(1, n)) .* (pbest_s(i, :) - subpop(i, :)) + abs(randn(1, n)) .* (lbest - subpop(i, :));

                        temp_s = subpop(i, :);

                        % Update the position of each particle
                        subpop(i, :) = subpop(i, :) + vel(i, :);

                        % Handle the elements of the particle which violate the boundary
                        tag1 = find(subpop(i, :)<min_var);
                        tag2 = find(subpop(i, :)>max_var);

                        if length(tag1)>0
                            subpop(i, tag1) = (temp_s(1, tag1) + min_var(1, tag1)) .* 0.5;
                        end

                        if length(tag2)>0
                            subpop(i, tag2) = (temp_s(1, tag2) + max_var(1, tag2)) .* 0.5;
                        end

                    end

                end

                p_              =  [p_; subpop];
                pbest_      =  [pbest_; pbest_s];
                fitPbest_   = [fitPbest_; fitPbest_s];

                k = k + 1;

            end

            p = [p; p_];
            pbest = [pbest; pbest_];
            fitPbest = [fitPbest; fitPbest_];

            % Evaluate the swarm
            fit = fitness(p, problem, delta, A);

            % Update the pbest of each particle based on the feasibility
            % criterion
            for i = 1:popsize

                if (fitPbest(i, 2) == 0 & fit(i, 2) == 0)

                    if fitPbest(i, 1)>fit(i, 1)
                        pbest(i, :) = p(i, :);
                        fitPbest(i, :) = fit(i, :);
                    end

                else

                    if fitPbest(i, 2)>fit(i, 2)
                        pbest(i, :) = p(i, :);
                        fitPbest(i, :) = fit(i, :);
                    end

                end

            end

            gen = gen + 1;

        end

        % Record the best results and the feasibility proportion
        feasiIndex = find(fitPbest(:, 2) == 0);
        if ~isempty(feasiIndex)
            bestResults = [bestResults min(fitPbest(feasiIndex, 1))];
            feasiPro = [feasiPro size(feasiIndex, 1)/popsize];
        else
            bestResults = [bestResults inf];
            feasiPro = [feasiPro 0];
        end

        time = time + 1;

    end

    % Show the best results and the feasibility proportion
    sort(bestResults)
    mean(bestResults)
    std(feasiPro)


toc;

% profview
