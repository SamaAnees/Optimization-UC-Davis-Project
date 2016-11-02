function rrange = return_range(r, Sig, num)
cvx_begin quiet
    R = chol(Sig);
    variable x(19)
    minimize(norm(R*x,2));

    subject to
    ones(1,19)*x == 1,
    x >= 0;
cvx_end

var = sqrt(x' * Sig * x);
exp = r' * x;

disp('Minimum risk:');
disp(x)

disp('Min expected return:')
disp(exp)

disp('Standard deviation min portfolio:')
disp(var)

cvx_begin quiet
	variable x(19)
    maximize(r' * x);
    subject to
    ones(1,19)*x == 1,
    x >= 0;
cvx_end
                    
maxU = cvx_optval;
                    
disp('Max return portfolio:')
disp(x);
rrange = linspace(exp, maxU, num);
                    
                    
end