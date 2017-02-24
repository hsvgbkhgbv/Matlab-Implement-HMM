clear;
N  = 100;         % number of sequences
T  = 100;        % length of the sequence
pi = [0.5; 0.5]; % inital probability pi_1 = 0.5 and pi_2 =0.5

%%two states hence A is a 2X2 matrix 
A  = [0.4 0.6 ; 0.4 0.6 ];         %p(y_t|y_{t-1})


%%one dimensional Gaussians 

E.mu    =[ .1   5 ]; %%the means of each of the Gaussians
E.sigma2=[ .4  .8 ]; %%the variances


[ Y, S ] = HmmGenerateData(N, T, pi, A, E, 'normal'); 

%%Y is the set of generated observations 
%%S is the set of ground truth sequence of latent vectors 

%Initilisation
B_e.mu = [0.5  7];
B_e.sigma2 = [0.6  2];

A_e = [ 0.6 0.4;
    0.2 0.8];

pi_e = [0.4;0.6];

%Set tolerance
tol = 0.0001;

%Assign value for comparison
A_ = A_e;
B_ = B_e;
pi_ = pi_e;

%Do EM steps to optimise parameters
for epoch=1:1000
    [E_z, E_z_z] = E_step(pi_e, A_e, B_e, Y, 'Continuous');
    [pi_e, A_e, B_e] = M_step(E_z, E_z_z, pi_e, A_e, B_e, Y, 'Continuous');
    if sum(sum(abs(A_ - A_e))) < tol && sum(sum(abs(B_.mu - B_e.mu))) < tol && sum(abs(pi_ - pi_e)) < tol && sum(sum(abs(B_.sigma2 - B_e.sigma2))) < tol
        break;
    else
        A_ = A_e;
        B_ = B_e;
        pi_ = pi_e;
    end
end

%Use vertibi to decode
[Y_e] = Viterbi(A_e, B_e, pi_e, Y, 'Continuous');

%Compute accuracy of results
accuracy = (N * T - sum(sum(abs(Y_e - S)))) / (N * T);