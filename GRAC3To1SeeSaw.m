% a function for generalized RAC, it takes as input a boolean array of length seven to mark
% which of the seven functions we need to maximize, the array specifies which
% of the seven functions is maxed [a0,a1,a2,a0+a1,a0+a2,a1+a2,a0+a1+a2]
function [vstepM,Rho,M] = GRAC3To1SeeSaw(a)
% preparing eight density matrices, indexed by three uniformly distributed
% bits
    Rho = cell(2,2,2); % for storing resultant states
    Rho_sdp = cell(2,2,2); % for storing sdp variables for states
    Fs = [];
    for a0 = 1:2
        for a1 = 1:2
            for a2 = 1:2
                Rho_sdp{a0,a1,a2} = sdpvar(2,2,'hermitian','complex'); 
                Fs =  [Fs; Rho_sdp{a0,a1,a2} >= 0; trace(Rho_sdp{a0,a1,a2}) == 1]; % positivity and unity trace
                Rho{a0,a1,a2} = RandomDensityMatrix(2);
            end
        end
    end

% preparing seven binary outcome measurements
% bits
    M = cell(7,2); % for storing resultant measurements
    M_sdp = cell(7,2); % for storing sdp variables for measurements
    Fm = [];
    for y = 1:7
        sum = 0;
        R = RandomPOVM(2,2); % random povm initialization for the first iteration
        for k = 1:2
            M_sdp{y,k} = sdpvar(2,2,'hermitian','complex');
            Fm = [Fm; M_sdp{y,k} >= 0;]; % positivity
            sum = sum + M_sdp{y,k};
            M{y,k} = R{k}; 
        end
        Fm = [Fm; sum == eye(2);]; % completeness
    end
    vstepM = 100; vstepS = 0; % declaring loop parameters 
    while (abs(vstepM-vstepS) >= 0.000000001)
        % sdp for states
        stepS = GRAC3To1Success(Rho_sdp,M,a); % the success metric
        % the state optimization step
        diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing states for the measurement sdp.
        for a0 = 1:2
            for a1 = 1:2
                for a2 = 1:2
                    Rho{a0,a1,a2} = value(Rho_sdp{a0,a1,a2});
                end
            end
        end
        % extracting the value of the objective function for the state step
        vstepS = value(stepS)
        % sdp for measurements.
        stepM = GRAC3To1Success(Rho,M_sdp,a); % the success metric
        % the measurement optimization step
        diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing states for the state sdp.
        for y = 1:7
            for k = 1:2
                M{y,k} = value(M_sdp{y,k});
            end
        end
         % extracting the value of the objective function for the measurement step
        vstepM = value(stepM)
    end