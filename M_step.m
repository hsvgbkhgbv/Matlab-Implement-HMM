%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function needs inputs
%pi = Kx1
%A = KxK
%B = Kx|obsevations|
%Y = NxT
%E_z = TxKxN
%E_z_z = TxKx(T-1)xN
%Model = 'Discrete' or 'Continuous'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function outputs
%pi_ = Kx1
%A_ = KxK
%B_ = Kx|obsevations|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Please input the correct model type!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pi_, A_, B_] = M_step(E_z, E_z_z, pi, A, B, Y, Model)
    
    %Initialisation
    N = size(Y, 1);
    T = size(Y, 2);
    K = length(pi);
    
    % Maximisation Step

    if strcmp(Model, 'Discrete')
        % Update emission matrix B
        B_ = zeros(size(B));
        
        for j=1:size(B, 2)
            for i=1:N
                for t=1:T
                    for k=1:K
                        if Y(i, t) == j     
                            B_(k, j) = B_(k, j) + E_z(t, k, i);  
                        end
                    end
                end
            end
        end
    
        normB = zeros(K, 1);
    
        for i=1:N
            normB = normB + reshape(sum(E_z(:, :, i), 1), [K, 1]);
        end

        normB = repmat(normB, 1, size(B_, 2));

        B_ = B_ ./ normB;

    elseif strcmp(Model, 'Continuous')
        % Update parameters of Gaussain Distribution
        B_.mu = [0  0];
        B_.sigma2 = [0  0];
    
        for t=1:T
            for i=1:N
                for k=1:K
                    B_.mu(k) = B_.mu(k) + E_z(t, k, i) * Y(i, t);
                    B_.sigma2(k) = B_.sigma2(k) + E_z(t, k, i) * ((Y(i, t) - B.mu(k)) * (Y(i, t) - B.mu(k))');
                end
            end
        end
    
        normB = zeros(K, 1);
    
        for i=1:N
            normB = normB + sum(E_z(:, :, i), 1)';
        end
        
        for k=1:K
            B_.mu(k) = B_.mu(k) / normB(k);
            B_.sigma2(k) = B_.sigma2(k) / normB(k);
        end
    end
    
    % Update pi
    pi_ = zeros(K, 1);

    for i=1:N
        pi_ = pi_ + E_z(1, :, i)';
    end

    normPi = 0;

    for i=1:N
        normPi = normPi + sum(E_z(1, :, i));
    end

    pi_ = pi_ / normPi;

    % Update transition matrix A
    A_ = zeros(size(A));

    for j=1:K
        for k=1:K
            for t=2:T
                for i=1:N
                    A_(j, k) = A_(j, k) + E_z_z(j, k, t - 1, i);
                end
            end
        end
    end

    normA = zeros(K, 1);

    for i=1:N
        for t=2:T
            normA = normA + sum(E_z_z(:, :, t - 1, i), 2);
        end
    end

    normA = repmat(normA, 1, K);
    
    A_ = A_ ./ normA;
end