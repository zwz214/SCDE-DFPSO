%**************************************************************************************************
%Author: Yong Wang
%Last Edited: July 1, 2010
%Email: ywang@csu.edu.cn
%Reference: Differential Evolution with Composite Trial Vector Generation Strategies
%                                                    and Control Parameters
%                           IEEE Transactions on Evolutionary Computation, Accepted
%**************************************************************************************************
function [position,value]= CoDE(benchmark_func,num_particle,dimension,runs,MaxDT,max_FES,func_num)
% clc;
% clear all;
% tic;

format long;
format compact;
y=(-1400:100:1400);
% 'CoDE'

% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.
% problemSet = [1 : 6 8 : 24];

%     problem = problemSet(problemIndex);
    problem = func_num;

    % Define the dimension of the problem
    n = dimension;
    lu = [-100 * ones(1, n); 100 * ones(1, n)];
%     switch problem
% 
%         case 1
% 
%             % lu: define the upper and lower bounds of the variables
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             % Load the data for this test function
%             load sphere_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 2
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load schwefel_102_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 3
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load high_cond_elliptic_rot_data
%             A = []; a = []; alpha = []; b = [];
% 
%             if n == 2, load elliptic_M_D2,
%             elseif n == 10, load elliptic_M_D10,
%             elseif n == 30, load elliptic_M_D30,
%             elseif n == 50, load elliptic_M_D50,
%             end
% 
%         case 4
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load schwefel_102_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 5
% 
%             lu = [-100 * ones(1,n); 100 * ones(1, n)];
%             load schwefel_206_data
%             M = []; a = []; alpha = []; b = [];
% 
%         case 6
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load rosenbrock_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 7
% 
%             lu = [0 * ones(1, n); 600 * ones(1, n)];
%             load griewank_func_data
%             A = []; a = []; alpha = []; b = [];
% 
%             c = 3;
%             if n == 2, load griewank_M_D2,
%             elseif n == 10, load griewank_M_D10,
%             elseif n == 30, load griewank_M_D30,
%             elseif n == 50, load griewank_M_D50,
%             end
% 
%         case 8
% 
%             lu = [-32 * ones(1, n); 32 * ones(1, n)];
%             load ackley_func_data
%             A = []; a = []; alpha = []; b = [];
% 
%             if n == 2, load ackley_M_D2,
%             elseif n == 10, load ackley_M_D10,
%             elseif n == 30, load ackley_M_D30,
%             elseif n == 50, load ackley_M_D50,
%             end
% 
%         case 9
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load rastrigin_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 10
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load rastrigin_func_data
%             A = []; a = []; alpha = []; b = [];
%             if n == 2, load rastrigin_M_D2,
%             elseif n == 10, load rastrigin_M_D10,
%             elseif n == 30, load rastrigin_M_D30,
%             elseif n == 50, load rastrigin_M_D50,
%             end
% 
%         case 11
% 
%             lu = [-0.5 * ones(1, n); 0.5 * ones(1, n)];
%             load weierstrass_data
%             A = []; a = []; alpha = []; b = [];
%             if n == 2, load weierstrass_M_D2,,
%             elseif n == 10, load weierstrass_M_D10,
%             elseif n == 30, load weierstrass_M_D30,
%             elseif n == 50, load weierstrass_M_D50,
%             end
% 
%         case 12
% 
%             lu = [-pi * ones(1, n); pi * ones(1, n)];
%             load schwefel_213_data
%             A = []; M = []; o = [];
% 
%         case 13
% 
%             lu = [-3 * ones(1, n); 1 * ones(1, n)];
%             load EF8F2_func_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 14
% 
%             lu = [-100 * ones(1, n); 100 * ones(1, n)];
%             load E_ScafferF6_func_data
%             if n == 2, load E_ScafferF6_M_D2,,
%             elseif n == 10, load E_ScafferF6_M_D10,
%             elseif n == 30, load E_ScafferF6_M_D30,
%             elseif n == 50, load E_ScafferF6_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 15
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func1_data
%             A = []; M = []; a = []; alpha = []; b = [];
% 
%         case 16
% 
%             lu = [-5 * ones(1,n); 5 * ones(1, n)];
%             load hybrid_func1_data
%             if n == 2, load hybrid_func1_M_D2,
%             elseif n == 10, load hybrid_func1_M_D10,
%             elseif n == 30, load hybrid_func1_M_D30,
%             elseif n == 50, load hybrid_func1_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 17
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func1_data
%             if n == 2, load hybrid_func1_M_D2,
%             elseif n == 10, load hybrid_func1_M_D10,
%             elseif n == 30, load hybrid_func1_M_D30,
%             elseif n == 50, load hybrid_func1_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 18
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func2_data
%             if n == 2, load hybrid_func2_M_D2,
%             elseif n == 10, load hybrid_func2_M_D10,
%             elseif n == 30, load hybrid_func2_M_D30,
%             elseif n == 50, load hybrid_func2_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 19
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func2_data
%             if n == 2, load hybrid_func2_M_D2,
%             elseif n == 10, load hybrid_func2_M_D10,
%             elseif n == 30, load hybrid_func2_M_D30,
%             elseif n == 50, load hybrid_func2_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 20
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func2_data
%             if n == 2, load hybrid_func2_M_D2,
%             elseif n == 10, load hybrid_func2_M_D10,
%             elseif n == 30, load hybrid_func2_M_D30,
%             elseif n == 50, load hybrid_func2_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 21
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func3_data
%             if n == 2, load hybrid_func3_M_D2,
%             elseif n == 10, load hybrid_func3_M_D10,
%             elseif n == 30, load hybrid_func3_M_D30,
%             elseif n == 50, load hybrid_func3_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 22
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func3_data
%             if n == 2, load hybrid_func3_HM_D2,
%             elseif n == 10, load hybrid_func3_HM_D10,
%             elseif n == 30, load hybrid_func3_HM_D30,
%             elseif n == 50, load hybrid_func3_HM_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 23
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func3_data
%             if n == 2, load hybrid_func3_M_D2,
%             elseif n == 10, load hybrid_func3_M_D10,
%             elseif n == 30, load hybrid_func3_M_D30,
%             elseif n == 50, load hybrid_func3_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 24
% 
%             lu = [-5 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func4_data
%             if n == 2, load hybrid_func4_M_D2,
%             elseif n == 10, load hybrid_func4_M_D10,
%             elseif n == 30, load hybrid_func4_M_D30,
%             elseif n == 50, load hybrid_func4_M_D50,
%             end
%             A = []; a = []; alpha = []; b = [];
% 
%         case 25
% 
%             lu = [2 * ones(1, n); 5 * ones(1, n)];
%             load hybrid_func4_data
%             if n == 2, load hybrid_func4_M_D2,
%             elseif n == 10, load hybrid_func4_M_D10,
%             elseif n == 30, load hybrid_func4_M_D30,
%             elseif n == 50, load hybrid_func4_M_D50,
%             end
%            A = []; a = []; alpha = []; b = [];
% 
%     end

    % Record the best results
    outcome = [];
    % Main body
    popsize = num_particle;

    time = 1;
    fitold=1;
    % The total number of runs
    totalTime = runs;
    Pg=zeros(1,max_FES);    %迭代过程中的最优适应度值 
    while time <= totalTime
  
        rand('seed', sum(100*clock));

        % Initialize the main population
        p = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

        % Evaluate the objective function values
        fit = benchmark_func(p', problem);
        fit=fit';
        % Record the number of function evaluations (FES)
        FES = popsize;
        m=0;
        while FES < max_FES
            m=m+1;
            pTemp = p;
            fitTemp = fit;
            
            % uSet: the set of trial vectors
            uSet = zeros(3 * popsize, n);

            for i = 1 : popsize

                % the three control parameter settings
                F    = [1.0 1.0 0.8];
                CR = [0.1 0.9 0.2];

                % Uniformly and randomly select one of the control
                % parameter settings for each trial vector generation strategy
                paraIndex = floor(rand(1, 3) * length(F)) + 1;

                % Generate the trail vectors
                u = generator(p, lu, i, F, CR, popsize, n, paraIndex);

                uSet(i * 3 - 2 : 3 * i, :) = u;

                FES = FES + 3;

            end

            % Evaluate the trial vectors
            fitSet = benchmark_func(uSet', problem);
            fitSet=fitSet';
            for i = 1 : popsize

                % Choose the best trial vector from the three trial vectors
                [~, minID] = min(fitSet(3 * i - 2 : 3 * i, :));
                bestInd = uSet(3 * (i - 1) + minID, :);
                bestIndFit = fitSet(3 * (i - 1) + minID, :);
                % Choose the better one between the trial vector and the
                % target vector
                if fit(i) >= bestIndFit
                    pTemp(i, :) = bestInd;
                    fitTemp(i, :) = bestIndFit;

                end
            end
            p = pTemp;
            fit = fitTemp;
            if time==1&&FES<=max_FES
                Pg(1,fitold:FES)=min(fit)-y(func_num);
            end
            fitold=FES;
         if FES>=max_FES
            break;
         end
        end
        outcome = [outcome min(fit)-y(func_num)];       
        time = time + 1;
    end
   position=Pg;
   value=outcome;
   
% toc;
