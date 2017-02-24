%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function needs inputs
%pi = Kx1
%A = KxK
%B = Kx|obsevations|
%S = 1xT
%C = Tx1
%Model = 'Discrete' or 'Continuous'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function outputs
%beta_ = TxK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Please input the correct model type!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [beta_] = beta(pi, A, B, S, C, Model)

    % initialise beta_
    K = length(pi);
    T = length(S);
    beta_ = zeros(T, K);
    
    % Compute \beta(z_{K})
    beta_(T, :) = ones(1, K);
    
    if strcmp(Model, 'Discrete')
        % Compute \beta(z_{1:(K - 1)})
        for t=(T - 1):-1:1
            beta_(t, :) = (A * (beta_(t + 1, :)' .* B(:, S(t + 1))) / C(t + 1))';
        end
    elseif strcmp(Model, 'Continuous')
        % Compute \beta(z_{1:(K - 1)})
        for t=T - 1:-1:1
            beta_(t, :) = (A * (beta_(t + 1, :)' .* [normpdf(S(t + 1), B.mu(1), sqrt(B.sigma2(1))); normpdf(S(t + 1), B.mu(2), sqrt(B.sigma2(2)))]) / C(t + 1))';
        end
    end
end