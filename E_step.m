%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function needs inputs
%pi = Kx1
%A = KxK
%B = Kx|obsevations|
%Y = NxT
%Model = 'Discrete' or 'Continuous'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function outputs
%E_z = TxKxN
%E_z_z = TxKx(T-1)xN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Please input the correct model type!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_z, E_z_z] = E_step(pi, A, B, Y, Model)

    %Expectation Step

    % Initialise the expectation matrix
    K = length(pi);
    N = size(Y, 1);
    T = size(Y, 2);
    E_z = zeros(T, K, N);
    E_z_z = zeros(K, K, T - 1, N);
    
    
    for i=1:N
        %Get Filtering result p(z_t | x_{1:t})
        [alpha_, C] = alpha(pi, A, B, Y(i, :), Model);
        
        %Get \beta(t)
        beta_ = beta(pi, A, B, Y(i, :), C, Model);
        
        % Compute Smoothing result p(z_t | x_1, x_2, ..., x_T)
        E_z(:, :, i) = alpha_ .* beta_;
        
        % Compute p(z_t, z_{t - 1} | x_1, x_2, ..., x_T)
        for t=2:T
            alphaTarget = alpha_(t - 1, :)';
            betaTarget = beta_(t, :);
            if strcmp(Model, 'Discrete')
                b_ = B(:, Y(i, t))';
            elseif strcmp(Model, 'Continuous')
                b_ = [normpdf(Y(i, t), B.mu(1), sqrt(B.sigma2(1))), normpdf(Y(i, t), B.mu(2), sqrt(B.sigma2(2)))];
            end
            E_z_z(:, :, t - 1, i) = A .* (alphaTarget * (b_ .* betaTarget)) / C(t);
        end
    end
end