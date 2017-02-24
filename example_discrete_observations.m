clear;
N  = 100;         % number of sequences
T  = 100;        % length of the sequence
pi = [0.5; 0.5]; % inital probability pi_1 = 0.5 and pi_2 =0.5

%%two states hence A is a 2X2 matrix 
A  = [0.4 0.6 ; 0.4 0.6 ];         %p(y_t|y_{t-1})

%%alphabet of 6 letters (e.g., a die with 6 sides) E(i,j) is the
E = [1/6 1/6 1/6 1/6 1/6 1/6;      %p(x_t|y_{t}) 
    1/10 1/10 1/10 1/10 1/10 1/2];

[ Y, S ] = HmmGenerateData(N, T, pi, A, E ); 

%%Y is the set of generated observations 
%%S is the set of ground truth sequence of latent vectors 


%Initialisation
B_e = [0.17 0.14 0.13 0.16 0.2 0.2;
    0.12 0.08, 0.16 0.04 0.2 0.4];

pi_e = [0.4; 0.6];

A_e = [0.3 0.7;
    0.5 0.5];

%Set tolerance
tol = 0.0001;

%Assign values for comparison
A_ = A_e;
B_ = B_e;
pi_ = pi_e;

%Execute EM step to optimise parameters
for epoch=1:100
    [E_z, E_z_z] = E_step(pi_e, A_e, B_e, Y, 'Discrete');
    [pi_e, A_e, B_e] = M_step(E_z, E_z_z, pi_e, A_e, B_e, Y, 'Discrete');
    if sum(sum(abs(A_ - A_e))) < tol && sum(sum(abs(B_ - B_e))) < tol && sum(abs(pi_ - pi_e)) < tol
        break;
    else
        A_ = A_e;
        B_ = B_e;
        pi_ = pi_e;
    end
end

%Use Viterbi to decode
[Y_e] = Viterbi(A_e, B_e, pi, Y, 'Discrete');

%Compute accuracy of results
accuracy = (N * T - sum(sum(abs(Y_e - S)))) / (N * T);