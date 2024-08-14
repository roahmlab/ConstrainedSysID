%% Numerical Regrouping using a QR decomposition

% Algorithm based on :
% Gautier et al. "Numerical calculation of the base inertial parameters of robots", 1990

% the inputs are the Observation Matrix which consists of multiple Regressor
% Matrices Y : [Y1, Y2, ..., Ynt]^T with 1...nt measurements
% and the Inertial Parameters InParam as symbolic variables.

% Inputs:
% ObservationMatrix: dimension: nq*nt x p , type : double

% Outputs:
% Selection-matrices Aid and Ad (can be used by multiplication on the Observation matrix and Inertial parameter vector, to get the independent (Aid) or dependent (Ad) parts)
% beta are the regrouped parameters or base inertial parameters
% Inverse mapping Ginv (InParam = Ginv*[beta; zeros(dim_d,1)])
% dim_id and dim_d are dimension of the independent and dependent parameter vectors
% Kd is a matrix used for computing Ginv and beta.

% Note:
% usually the dependent inertial parameters are set to zero because of
% simplicity and to get a solution for the regrouped parameters
% beta.

% Author: Daniel Haugk, 2024, University of Michigan

function [Aid, Ad, Kd, beta, Ginv, dim_id, dim_d] = QRDecomposition(ObservationMatrix,InParam)

% assign Observation matrix to variable
W = ObservationMatrix;

% get number of parameters
[~, p] = size(W);

% get permutation vector M which contains the index of each independant column (Parameter) as its first "rankW" elements
[~, ~, M] = qr(W,0);

% rank estimation
rankW = rank(W); 

% get index of each independant parameter and sort them
idx=sort(M(1:rankW));

% get each index of dependant columns (Parameters), by checking which elements are in the independant index set
j = 1;
for i = 1:p

    if ~ismember(i,idx)
        idx_(j) = i;
        j = j + 1;
    end

end

% sort dependant index set
idx_ = sort(idx_);

% intialize selection matrices
Aid = zeros(p,rankW);
Ad = zeros(p,p-rankW);

% get selection matrix for indepenant parameters
for i = 1:rankW

    Aid(idx(i),i) = 1; 

end

% get selection matrix for dependant parameters
for i = 1:p-rankW

    Ad(idx_(i),i) = 1;

end

% combine both selection matrix to get full selection matrix
A = [Aid Ad];

% get QR-decomposition matrix R for computation of Kd (see Gautier paper) 
[~ , R_] = qr(W*A);

% split R matrix into an indepant and dependant part
R1 = R_(1:rankW,1:rankW);
R2 = R_(1:rankW,rankW+1:end);

% compute Kd
Kd = R1\R2;

% get rid of extremely small numbers for more numerical stability
Kd(abs(Kd)<sqrt(eps)) = 0;

% get independant parameters
pi_id = Aid.'*InParam;

% get dependant parameters
pi_d = Ad.'*InParam;

% compute matrix for G/Ginv
KG_ = [eye(length(pi_id)) -Kd; zeros(length(pi_d),length(pi_id)) eye(length(pi_d))];

% compute base inertial parameters
beta = pi_id + Kd*pi_d;

% compute Ginv
Ginv = A*KG_;

% dimension of independant parameters
dim_id = length(pi_id);

% dimension of dependant parameters
dim_d = length(pi_d);
    
end