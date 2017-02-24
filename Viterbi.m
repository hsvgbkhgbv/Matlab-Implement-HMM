%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function needs inputs
%pi = Kx1
%A = KxK
%B = Kx|obsevations|
%Y = NxT
%Model = 'Discrete' or 'Continuous'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function outputs
%Y_e = NxT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Please input the correct model type!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y_e] = Viterbi(A, B, pi, Y, Model)

    %Viterbi Algorithm initialisation
    K = size(A, 1);
    T = size(Y, 2);
    N = size(Y, 1);
    Y_e = zeros(N, T);
   
    for i=1:N
        scores = zeros(K, T);
        latent_vars = zeros(K, T - 1);
        if strcmp(Model, 'Discrete')
            initial = log2(pi) + log2(B(:, Y(i, 1)));
            scores(:, 1) = initial;
            for t=2:T
                iter1 = log2(A(:, 1)) + [log2(B(1, Y(i, t)));
                    log2(B(1, Y(i, t)))] + scores(:, t - 1);
                iter2 = log2(A(:, 2)) + [log2(B(2, Y(i, t)));
                    log2(B(2, Y(i, t)))] + scores(:, t - 1);
                [M1, I1] = max(iter1);
                [M2, I2] = max(iter2);
                latent_vars(:, t - 1) = [I1; I2];
                scores(:, t) = [M1; M2];
            end
        elseif strcmp(Model, 'Continuous')
            initial = log2(pi) + [log2(normpdf(Y(i, 1), B.mu(1), sqrt(B.sigma2(1))));
                log2(normpdf(Y(i, 1), B.mu(2), sqrt(B.sigma2(2))))];
            scores(:, 1) = initial;
            for t=2:T
                iter1 = log2(A(:, 1)) + [log2(normpdf(Y(i, t), B.mu(1), sqrt(B.sigma2(1))));
                log2(normpdf(Y(i, t), B.mu(1), sqrt(B.sigma2(1))))] + scores(:, t - 1);
                iter2 = log2(A(:, 2)) + [log2(normpdf(Y(i, t), B.mu(2), sqrt(B.sigma2(2))));
                log2(normpdf(Y(i, t), B.mu(2), sqrt(B.sigma2(2))))] + scores(:, t - 1);
                [M1, I1] = max(iter1);
                [M2, I2] = max(iter2);
                latent_vars(:, t - 1) = [I1; I2];
                scores(:, t) = [M1; M2];
            end
        end
        [~, Y_e(i, T)] = max(scores(:, T));
        for t=T-1:-1:1
            Y_e(i, t) = latent_vars(Y_e(i, t + 1), t);
        end
    end
end