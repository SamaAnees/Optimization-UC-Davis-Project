function [Y, rates, sigs] = efficient_frontier(r, Sig, num)

n = length(r);
[rrange] = return_range(r, Sig, num);

Y = zeros(n, num);
for i = 1:num
    cvx_begin quiet
        R = chol(Sig);
        variable x(19)
        minimize(norm(R*x, 2));
        subject to
        ones(1,19)*x == 1,
        r'*x >= rrange(i)
        x >= 0;
    cvx_end
    Y(:,i) = x;
    rates(i) = r'*x;
    sigs(i) = sqrt(x'*Sig*x);
end

rates = rates';
sigs = sigs';

end