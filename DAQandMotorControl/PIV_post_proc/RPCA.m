%% Robust Principle Component Analysis
% From Brunton and Kutz, Data-driven Science and Engineering, 2nd ed., 2022
% Updated 20220123

function [L,S] = RPCA(X)
    
    tStart = tic;
    [n1,n2] = size(X);
    mu = n1*n2/(4*sum(abs(X(:))));
    lambda = 1/sqrt(max(n1,n2));
    thresh = 1e-7*norm(X,'fro');
    
    L = zeros(size(X));
    S = zeros(size(X));
    Y = zeros(size(X));
    count = 0;
    
    while((norm(X-L-S,'fro')>thresh) && (count<1000))
        L = SVT(X-S+(1/mu)*Y,1/mu);
        S = shrink(X-L+(1/mu)*Y,lambda/mu);
        Y = Y + mu*(X-L-S);
        count = count + 1;
    end
    tEnd = toc(tStart);
    disp(['RPCA calculation took ', num2str(tEnd/60/60), ' hours.'])
    
    % Nested functions
    
    function out = SVT(X, tau)
        [U,S,V] = svd(X,'econ');
        out = U*shrink(S,tau)*V';
    end
    
    function out = shrink(X, tau)
        out = sign(X).*max(abs(X) - tau, 0);
    end

end