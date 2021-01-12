% a function for generalized RAC, it takes as input a boolean array of length seven to mark
% which of the seven functions we need to maximize, the array specifies which
% of the seven functions is maxed [a0,a1,a2,a0+a1,a0+a2,a1+a2,a0+a1+a2]
function [vstepMB,Rho,MA,MB] = GRAC3To1EACCSeeSaw(a,q)
% preparing the potentially entangled state 
    Rho_sdp = sdpvar(q^2,q^2,'hermitian','complex'); % for storing sdp variables for states
    Fs = [Rho_sdp >= 0; trace(Rho_sdp) == 1];
% preparing Alice's eight binary outcome measurements 
    MA = cell(2,2,2,2); % for storing resultant Alice's measurements
    MA_sdp = cell(2,2,2,2); % for storing sdp variables for Alice's measurements
    Fma = [];
    for a0 = 1:2
        for a1 = 1:2
            for a2 = 1:2
                sum=0;
                R = RandomPOVM(q,2);
                for m = 1:2
                    MA_sdp{a0,a1,a2,m} = sdpvar(q,q,'hermitian','complex'); 
                    sum = sum+MA_sdp{a0,a1,a2,m};
                    Fma =  [Fma; MA_sdp{a0,a1,a2,m} >= 0; ]; % positivity
                    MA{a0,a1,a2,m} = R{m}; 
                end
                Fma =  [Fma; sum == eye(q); ]; % completeness
            end
        end
    end
% preparing Bob's fourteen binary outcome measurements
% bits
    MB = cell(7,2,2); % for storing resultant measurements
    MB_sdp = cell(7,2,2); % for storing sdp variables for measurements
    Fmb = [];
    for y = 1:7
        for m = 1:2
            sum = 0;
            R = RandomPOVM(q,2); % random povm initialization for the first iteration
            for k = 1:2
                MB_sdp{y,m,k} = sdpvar(q,q,'hermitian','complex');
                Fmb = [Fmb; MB_sdp{y,m,k} >= 0;]; % positivity
                sum = sum + MB_sdp{y,m,k};
                MB{y,m,k} = R{k}; 
            end
            Fmb = [Fmb; sum == eye(q);]; % completeness
        end
    end
        
        
        
    vstepMB = 7; vstepMB = 77; vstepS = 777; tol = 0.000000001; % declaring loop parameters 
    while (abs(vstepMB-vstepS) >= tol) ||  (abs(vstepMA-vstepS) >= tol) ||  (abs(vstepMB-vstepMA) >= tol)    
        % sdp for states
        stepS = GRAC3To1EACCSuccess(Rho_sdp,MA,MB,a); % the success metric
        % the state optimization step
        diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing state for the measurement sdp.
        Rho = value(Rho_sdp);
        vstepS = value(stepS)
        
        % sdp for Alice's measurements
        stepMA = GRAC3To1EACCSuccess(Rho,MA_sdp,MB,a); % the success metric
        % the state optimization step
        diagnostics = optimize([Fma;], -stepMA, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing state for the measurement sdp.
        for a0 = 1:2
            for a1 = 1:2
                for a2 = 1:2
                    for m = 1:2
                        MA{a0,a1,a2,m} = value(MA_sdp{a0,a1,a2,m});
                    end
                end
            end
        end
        % extracting the value of the objective function for the state step
        vstepMA = value(stepMA)
        
        % sdp for Bob's measurements.
        stepMB = GRAC3To1EACCSuccess(Rho,MA,MB_sdp,a); % the success metric
        % the measurement optimization step
        diagnostics = optimize([Fmb;], -stepMB, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing states for the state sdp.
        for y = 1:7
            for m = 1:2
                for k = 1:2
                    MB{y,m,k} = value(MB_sdp{y,m,k});
                end
            end
        end
         % extracting the value of the objective function for the measurement step
        vstepMB = value(stepMB)
    end