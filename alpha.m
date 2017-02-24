%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function needs inputs
%pi = Kx1
%A = KxK
%B = Kx|obsevations|
%S = 1xT
%Model = 'Discrete' or 'Continuous'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function outputs
%alpha_ = TxK
%C = Tx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Please input the correct model type!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha_, C] = alpha(pi, A, B, S, Model)

    % initialise alpha_
    K = size(pi, 1);
    T = length(S);
    alpha_ = zeros(T, K);
    C = zeros(T, 1);
    
    if strcmp(Model, 'Discrete')
        % Compute \alpha(z_1)
        
        C(1) = sum(pi .* B(:, 1));
        alpha_(1, :) = (pi .* B(:, 1) / C(1))';
        
        % Compute \alpha(z_{2:T})
        for t=2:T
            C(t) = sum(B(:, S(t)) .* (A' * alpha_(t - 1, :)'));
            alpha_(t, :) = (B(:, S(t)) .* (A' * alpha_(t - 1, :)') / C(t))';
        end

    elseif strcmp(Model, 'Continuous')
        % Compute \alpha(z_1)
        C(1) = sum(pi .* [normpdf(S(1), B.mu(1), sqrt(B.sigma2(1))); normpdf(S(1), B.mu(2), sqrt(B.sigma2(2)))]);
        alpha_(1, :) = (pi .* [normpdf(S(1), B.mu(1), sqrt(B.sigma2(1))); normpdf(S(1), B.mu(2), sqrt(B.sigma2(2)))] / C(1))';
    
        % Compute \alpha(z_{2:T})
        for t=2:T
            C(t) = sum([normpdf(S(t), B.mu(1), sqrt(B.sigma2(1))); normpdf(S(t), B.mu(2), sqrt(B.sigma2(2)))] .* (A' * alpha_(t - 1, :)'));
            alpha_(t, :) = ([normpdf(S(t), B.mu(1), sqrt(B.sigma2(1))); normpdf(S(t), B.mu(2), sqrt(B.sigma2(2)))] .* (A' * alpha_(t - 1, :)') / C(t))';
        end
    end
end