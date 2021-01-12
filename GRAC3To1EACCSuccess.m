function S = GRAC3To1EACCSuccess(Rho,MA,MB,a)
% for calculating the success metric given states and measurements, and the
% vector a specifying the questions of Bob 
    S=0;
    for a0 = 1:2
        for a1 = 1:2
            for a2 = 1:2
                for m = 1:2
                    if a(1) == 1 % a0
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},MB{1,m,a0})); 
                    end
                    if a(2) == 1 % a1
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},MB{2,m,a1}));
                    end
                    if a(3) == 1 % a2
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},MB{3,m,a2}));
                    end
                    if a(4) == 1 % a0+a1
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},MB{4,m,xor(a0-1,a1-1)+1}));
                    end
                    if a(5) == 1 % a0+a2
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},MB{5,m,xor(a0-1,a2-1)+1}));
                    end
                    if a(6) == 1 % a1+a2
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},M{6,m,xor(a1-1,a2-1)+1}));
                    end
                    if a(7) == 1 % a0+a2
                        S = S + trace(Rho*kron(MA{a0,a1,a2,m},MB{7,m,xor(xor(a1-1,a2-1),a0-1)+1}));
                    end
                end
            end
        end
    end
    S = real(S/(8*sum(a)));