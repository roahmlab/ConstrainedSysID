%% IRLS Algorithm for only inertial parameters without any normalization and weighting

% Inputs:

% n = number of joints, double
% na = number of actuated joints, double
% nu = number of unactuated joints, double
% na/nu_idx = vector of indices of actuated/unactuated joints, na/nu vector
% m = number of data points, double
% b = number of base inertial parameters, double
% d = number of dependant parameters, double
% W_ip = Observation matrix for inertial parameters, n*m x 10*n matrix 
% T = torque vector over m data points, n*m vector
% data = matrix containing pos., vel. and accel. for all joints over m data points, n x 3*m matrix
% Ginv = inverse of bijective map to get inertial parameters from base inertial parameters, 10*n x 10*n matrix
% Aid = selection matrix for indepenent parameters, 10*n x b matrix
% Ad = selection matrix for dependant and unidentifiable parameters, 10*n x (10*n-b) matrix
% Kd = regrouping transformation matrix, b x (10*n-b) matrix
% J = constraint jacobian function, nc x n matrix (nc = number of constraints)
% lb = lower bound on paramaters, p_full vector
% ub = upper bound on parameters p_full vector
% X_fm = motor, friction and offset parameters and friction exponent values, (p_full-9*n) vector
% X0_1 = initial condition for parameters, p_full vector

% Algorithm options:

% MS1 = number of multi searchs for first optimization, double
% regroup = include regouping or not, bool
% includeFMDynamics = include motor and friction dynamics, bool
% includeOffset = include offset as parameter, bool
% includeConstraints = systen has constraints, bool
% constraintVariant = choose method of constraint projection matrix, double (1 or 2)

% Output:

% X = vector of all parameters from optimization, p_full vector 
% Wfull = observation matrix for inertial, motor, friction and offset parameters, n*m x p_full matrix
% alphanew = new friction exponents from optimization, n vector

% Author: Daniel Haugk, 2024, University of Michigan

function [X, Wfull, alphanew] = SysIDAlgInParamOnlyNoWeights(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,X_fm,T,data,dataFull,W_ip,Ginv,Aid,Ad,Kd,X0_1,options)

% initialize variables
n = variables{1};
m = variables{2};
b = variables{3};
d = variables{4};

% initalize algorithm options
MS1 = AlgOptions{3};
regroup = AlgOptions{5};
includeOffset = AlgOptions{9};
includeFMDynamics = AlgOptions{10};
includeConstraints = AlgOptions{11};
constraintVariant = AlgOptions{12};

% compute number of all parameters with regrouping
if includeOffset
    b_full = b + 4*n;
else
    b_full = b + 3*n;
end

% check if constraints are included
if includeConstraints

    % get number of actuated links
    na = length(na_idx);

    if constraintVariant == 1
    % general case for all kinds of constraints

        for i = 1:m

            % compute constraint projection matrix
            J_temp = J(dataFull(:,i));
            K = eye(n) - J_temp(:,sort([na_idx nu_idx])).'*pinv(J_temp(:,sort([na_idx nu_idx])).');

            % compute constrained observation matrix and joint torques
            T((i-1)*n+1:i*n) = K*T((i-1)*n+1:i*n);
            W_ip((i-1)*n+1:i*n,:) = K*W_ip((i-1)*n+1:i*n,:);
        end

    elseif constraintVariant == 2
    % special case for fully actuated systems

        W_ip_temp = [];

        for i = 1:m

            % compute constraint projection matrix
            J_temp = J(dataFull(:,i));
            P = [eye(na); -J_temp(:,nu_idx)\J_temp(:,na_idx)];

            % compute constrained observation matrix on actuated joint level
            W_ip_ = P.'*W_ip((i-1)*n+1:i*n,:);
            W_ip_temp = [W_ip_temp;W_ip_];
        end

        % assign temporary observation matrix variable to observation matrix
        W_ip = W_ip_temp;

    end

else
    
    % set variables for constrained system to empty
    na = [];
    na_idx = [];
    nu_idx = [];
    J = [];

end

% compute T_hat if friction and motor dynamics are included and their parameters are provided
if includeFMDynamics
    % get motor, friction and offset observation matrix
    [Wfm,~,~] = Wextension(X_fm(1:n),n,na,na_idx,nu_idx,J,m,data,dataFull,includeOffset,includeConstraints,constraintVariant);
    
    % compute T_hat
    T_hat = T - Wfm*X_fm(n+1:end);
else
    % assign T_hat
    T_hat = T;
end

% check if regrouping is wanted
if regroup
    % get regrouped observation matrix
    Wb = W_ip*Aid;

    % compute regrouped bounds for optimization problem
    lowerB1 = Aid.'*lb(1:10*n) + Kd*Ad.'*lb(1:10*n);
    upperB1 = Aid.'*ub(1:10*n) + Kd*Ad.'*ub(1:10*n);

    % check if regrouped bounds have to be switched
    for j = 1:length(lowerB1)
        if lowerB1(j)>=upperB1(j)
            temp_lb = lowerB1(j);
            temp_ub = upperB1(j);
            lowerB1(j) = temp_ub;
            upperB1(j) = temp_lb;
        end
    end

    % regroup the initial conditions
    X0_ip = Aid.'*X0_1(1:10*n) + Kd*Ad.'*X0_1(1:10*n);

    % assign regrouped iniertial parameters
    X0_1 = X0_ip;   

else
    % get observation matrix
    Wb = W_ip;

    % get bounds for optimization problem
    lowerB1 = lb;
    upperB1 = ub;
end

% run optimization algorithm
[X, W_ip] = AlgInParamOnly(Wb,T_hat,n,b,b_full,d,MS1,lowerB1,upperB1,X0_1,regroup,Ginv,options);

% put together full observation matrix
Wfull = [W_ip Wfm];

% assign friction exponents
alphanew = X_fm(1:n);

end


%% Algorithm functions

% Algorithm without optimimzation over friction exponent
function [X, W_ip] = AlgInParamOnly(Wb,T_hat,n,b,b_full,d,MS1,lowerB1,upperB1,X0_1,regroup,Ginv,options)
   
     % assign observation matrix for inertial parameters
     W_ip = Wb;

     % initialize least squares function
     fun1 = @(X1) optimFun(X1, W_ip, T_hat);

     % formulate least squares problem for all parameters
     problem1 = createOptimProblem('fmincon','x0',X0_1,'objective',fun1,'lb',lowerB1,'ub',upperB1,'nonlcon',@(X1) lmiconDet(X1,n,Ginv,b,b_full,d,regroup),'options',options);
                
     % choose multi start option to increase the search space
     gs = MultiStart;

     % start optimization
     X = run(gs,problem1,MS1);
                
end


%% Least squares functions

% least squares function for all parameters excluding friction exponent
function [fmin, gradf] = optimFun(X, W, T)

    % least squares function
    fmin = norm(W*X - T)^2;
    
    % gradient of least squares function
    if nargout > 1
        gradf = 2*W.'*(W*X - T);
    end
 
end


%% Helper functions

% Computation of motor and friction observation matrix
% Wm: is the motor observation matrix
% Wf: is the friction and offset observtation matrix
% Wfm: is both combined
function [Wfm, Wm, Wf] = Wextension(alpha,nLinks,na,na_idx,nu_idx,J,mData,data,dataFull,includeOffset,includeConstraints,constraintVariant)

    % check if offset is included and preallocate friction regressor and observation matrix
    if includeOffset
        Wf = zeros(nLinks*mData,3*nLinks);
        Yf = zeros(nLinks,3*nLinks);
    else
        Wf = zeros(nLinks*mData,2*nLinks);
        Yf = zeros(nLinks,2*nLinks); 
    end

    % preallocate motor observation matrix
    Wm = zeros(nLinks*mData,nLinks);
    % preallocate motor regressor matrix
    Ym = zeros(nLinks,nLinks);
    
    % compute regressor matrices over all data points
    for i = 1:mData

        % compute regressor matrices
        for j = 1:nLinks

            % assign accelerations for each joint in motor regressor matrix : Fm = motorInertia*qdd -> dFm/dmotorInertia = qdd 
            Ym(j,j) = data(j,i + 2*mData);

            % check if offset is included and assign dFf/d(Fc,Fv,beta) for each joint, with Fc: coloumb friction coefficient, Fv: vicious friction coefficient, beta: offset
            if includeOffset
                Yf(j,3*(j-1)+1:3*j) = [sign(data(j,i + mData)) sign(data(j,i + mData))*abs(data(j,i + mData))^alpha(j) 1];
            else
                Yf(j,2*(j-1)+1:2*j) = [sign(data(j,i + mData)) sign(data(j,i + mData))*abs(data(j,i + mData))^alpha(j)];
            end
        end

        % assign regressor matrices to observation matrix
        Wm((i-1)*nLinks+1:i*nLinks,:) = Ym;
        Wf((i-1)*nLinks+1:i*nLinks,:) = Yf;

    end

    % check if constraints are included
    if includeConstraints
        
        if constraintVariant == 1
        % general case for all kinds of constraints
        
            for i = 1:mData
        
                % compute constraint projection matrix
                J_temp = J(dataFull(:,i));
                K = eye(nLinks) - J_temp(:,sort([na_idx nu_idx])).'*pinv(J_temp(:,sort([na_idx nu_idx])).');
        
                % compute constrained observation matrix for motor and friction dynamics
                Wm((i-1)*nLinks+1:i*nLinks,:) = K*Wm((i-1)*nLinks+1:i*nLinks,:);
                Wf((i-1)*nLinks+1:i*nLinks,:) = K*Wf((i-1)*nLinks+1:i*nLinks,:);
            end
        
        elseif constraintVariant == 2
        % special case for fully actuated systems
        
            Wf_temp = [];
            Wm_temp = [];
        
            for i = 1:mData
        
                % compute constraint projection matrix
                J_temp = J(dataFull(:,i));
                P = [eye(na); -J_temp(:,nu_idx)\J_temp(:,na_idx)];
        
                % compute constrained observation matrix for motor and friction dynamics on actuated joint level
                Wm_ = P.'*Wm((i-1)*nLinks+1:i*nLinks,:);
                Wf_ = P.'*Wf((i-1)*nLinks+1:i*nLinks,:);
                Wm_temp = [Wm_temp;Wm_];
                Wf_temp = [Wf_temp;Wf_];
            end
        
            % assign temporary observation matrix variable to observation matrix
            Wm = Wm_temp;
            Wf = Wf_temp;
        
        end
   end
        
   % get full observation matrix for motor and friction parameters
   Wfm = [Wm Wf];

end 


% ------------------------------------------------------------------------


% nonlinear LMI constraint
function [LMIdet , nonlinEqCon, DC, DCeq] = lmiconDet(X,nLinks,Ginv,b,b_full,d,regroup)
    
    %check if regrouping was done
    if regroup
        % get full inertial parameter vector
        phi = Ginv*[X(1:b); zeros(d,1)];
        % derivative of phi w.r.t. x for LMI-gradient
        dphidX = Ginv*[eye(b) zeros(b,b_full-b); zeros(d,b_full)];
        % preallocate gradient of LMI-matrix determinant
        DC = zeros(nLinks,10*nLinks);
    else
        % get full inertial parameter vector
        phi = X;
        % preallocate gradient of LMI-matrix determinant
        DC = zeros(nLinks,nLinks*10+3*nLinks);
    end

    % preallocate LMI-matrix
    LMIdet = zeros(nLinks,1);
    
    % compute derivative of LMI-matrix w.r.t. inertial parameters : dphi/d(LMI-matrix)
    if nargout > 2

        % get derivative of each element in the intertial parameters phi
        dphi1dLMI = [-0.5 0 0 0; 0 0.5 0 0; 0 0 0.5 0; 0 0 0 0];
        dphi2dLMI = [0 -1 0 0; -1 0 0 0; 0 0 0 0; 0 0 0 0];
        dphi3dLMI = [0 0 -1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0];
        dphi4dLMI = [0.5 0 0 0; 0 -0.5 0 0; 0 0 0.5 0; 0 0 0 0];
        dphi5dLMI = [0 0 0 0; 0 0 -1 0; 0 -1 0 0; 0 0 0 0];
        dphi6dLMI = [0.5 0 0 0;0 0.5 0 0; 0 0 -0.5 0; 0 0 0 0];
        dphi7dLMI = [0 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 0];
        dphi8dLMI = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 1 0 0];
        dphi9dLMI = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 1 0];
        dphi10dLMI = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1];

        % plug them together to get the full derivative matrix
        dphidLMI = [dphi1dLMI dphi2dLMI dphi3dLMI dphi4dLMI dphi5dLMI ...
                    dphi6dLMI dphi7dLMI dphi8dLMI dphi9dLMI dphi10dLMI];

    end

    % compute LMI-matrix
    for i = 1:nLinks

        % get inertia tensor
        I = [phi(10*(i-1)+1) phi(10*(i-1)+2) phi(10*(i-1)+3);
             phi(10*(i-1)+2) phi(10*(i-1)+4) phi(10*(i-1)+5)
             phi(10*(i-1)+3) phi(10*(i-1)+5) phi(10*(i-1)+6)];
    
        % compute all 4 parts of the LMI-matrix
        lmi11 = (trace(I)/2)*eye(3) - I;
        lmi12 = [phi(10*(i-1)+7) phi(10*(i-1)+8) phi(10*(i-1)+9)]';
        lmi21 = lmi12';
        lmi22 = phi(10*(i-1)+10);
    
        % plug them together to get the LMI-matrix
        LMI = [lmi11 lmi12; lmi21 lmi22];

        % get determinant of LMI-matrix (negative sign for positive definite check)
        LMIdet(i,1) = -det(LMI);

        % compute gradient of determinant of LMI-matrix
        if nargout > 2
           
            % compute gradient for each parameter seperately
            for j = 1:10
                % check if LMI is singular, if it is use pseudo inverse
                if rank(LMI) < 4
                    DC(i,j+(i-1)*10) = LMIdet(i,1)*trace(pinv(LMI)*dphidLMI(:,4*(j-1)+1:j*4));  
                else
                    DC(i,j+(i-1)*10) = LMIdet(i,1)*trace(LMI\dphidLMI(:,4*(j-1)+1:j*4));
                end
            end

            % no equality constraints
            DCeq = [];
        end

    end
    
    % get full LMI-gradient (with regrouping we have the chain-rule: dLMI/dX = dLMI/dphi * dphi/dX)
    if regroup
        DC = (DC*dphidX).';
    else
        DC = DC.';
    end

    % no equality constraints
    nonlinEqCon = [];

end