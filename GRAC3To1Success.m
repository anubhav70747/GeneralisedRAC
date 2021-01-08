function S = GRAC3To1Success(Rho,M,a)
% for calculating the success metric given states and measurements, and the
% vector a specifying the questions of Bob 
    S=0;
    for a0 = 1:2
        for a1 = 1:2
            for a2 = 1:2
                if a(1) == 1 % a0
                    S = S + trace(Rho{a0,a1,a2}*M{1,a0}); 
                end
                if a(2) == 1 % a1
                    S = S + trace(Rho{a0,a1,a2}*M{2,a1});
                end
                if a(3) == 1 % a2
                    S = S + trace(Rho{a0,a1,a2}*M{3,a2});
                end
                if a(4) == 1 % a0+a1
                    S = S + trace(Rho{a0,a1,a2}*M{4,xor(a0-1,a1-1)+1});
                end
                if a(5) == 1 % a0+a2
                    S = S + trace(Rho{a0,a1,a2}*M{5,xor(a0-1,a2-1)+1});
                end
                if a(6) == 1 % a1+a2
                    S = S + trace(Rho{a0,a1,a2}*M{6,xor(a1-1,a2-1)+1});
                end
                if a(7) == 1 % a0+a2
                    S = S + trace(Rho{a0,a1,a2}*M{7,xor(xor(a1-1,a2-1),a0-1)+1});
                end
            end
        end
    end
    S = real(S/(8*sum(a)));