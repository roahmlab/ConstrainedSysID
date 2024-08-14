%% IRLS Algorithm for all parameters including motor, friction and offset parameters without any normalization and weighting

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
% lba = lower bound on friction exponents, n vector
% uba = upper bound on friction exponents, n vector
% X0_1 = initial condition for parameters, p_full vector

% Algorithm options:

% tol = stopping criterion for while loops, double << 1
% MS1 = number of multi searchs for first optimization, double
% MS2 = number of multi searchs for second optimization, double
% regroup = include regouping or not, bool
% SearchAlpha = include optimization over friction exponents, bool
% includeOffset = include offset as parameter, bool
% includeConstraints = systen has constraints, bool
% constraintVariant = choose method of constraint projection matrix, double (1 or 2)

% Output:

% X = vector of all parameters from optimization, p_full vector 
% Wfull = observation matrix for inertial, motor, friction and offset parameters, n*m x p_full matrix
% alphanew = new friction exponents from optimization, n vector

% Author: Daniel Haugk, 2024, University of Michigan

function [X, Wfull, alphanew] = SysIDAlgFullNoWeights(variables,AlgOptions,na_idx,nu_idx,J,lb,ub,lba,uba,T,data,dataFull,W_ip,Ginv,Aid,Ad,Kd,X0_1,options)

% initialize variables
n = variables{1};
m = variables{2};
b = variables{3};
d = variables{4};

% initalize algorithm options
tol = AlgOptions{2};
MS1 = AlgOptions{3};
MS2 = AlgOptions{4};
regroup = AlgOptions{5};
SearchAlpha = AlgOptions{8};
includeOffset = AlgOptions{9};
includeConstraints = AlgOptions{11};
constraintVariant = AlgOptions{12};

% compute number of inertial parameters
p_ip = 10*n;

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

% check if regrouping is wanted
if regroup

    % get regrouped observation matrix
    Wb = W_ip*Aid;

    % compute regrouped bounds
    lowerB_ = Aid.'*lb(1:10*n) + Kd*Ad.'*lb(1:10*n);
    upperB_ = Aid.'*ub(1:10*n) + Kd*Ad.'*ub(1:10*n);

    % check if regrouped bounds have to be switched
    for j = 1:length(lowerB_)
        if lowerB_(j)>=upperB_(j)
            temp_lb = lowerB_(j);
            temp_ub = upperB_(j);
            lowerB_(j) = temp_ub;
            upperB_(j) = temp_lb;
        end
    end

    % get bounds for first optimization problem
    lowerB1 = [lowerB_;lb(10*n+1:end)];
    upperB1 = [upperB_;ub(10*n+1:end)];

    % get bounds for second optimization problem
    lowerB2 = [lba;lb(11*n+1:end)];
    upperB2 = [uba;ub(11*n+1:end)];

    % regroup the initial conditions
    X0_ip = Aid.'*X0_1(1:10*n) + Kd*Ad.'*X0_1(1:10*n);
    X0_fm = X0_1(10*n+1:end);

    % put togehter regrouped initial condition
    X0_1 = [X0_ip;X0_fm];

else
    % get observation matrix
    Wb = W_ip;

    % get bounds for first optimization problem
    lowerB1 = lb;
    upperB1 = ub;

    % get bounds for second optimization problem
    lowerB2 = [lba;lb(11*n+1:end)];
    upperB2 = [uba;ub(11*n+1:end)];
end


% check if optimization over friction exponent is included
if SearchAlpha
    [X, Wfull, alphanew] = AlgAlpha(Wb,T,tol,n,na,na_idx,nu_idx,m,p_ip,b,b_full,d,MS1,MS2,data,dataFull,J,lowerB1,upperB1,lowerB2,upperB2,X0_1,regroup,Ginv,includeOffset,includeConstraints,constraintVariant,options);
else
    [X, Wfull] = AlgNoAlpha(Wb,T,n,na,na_idx,nu_idx,m,b,b_full,d,MS1,data,dataFull,J,lowerB1,upperB1,X0_1,regroup,Ginv,includeOffset,includeConstraints,constraintVariant,options);
    alphanew = ones(n,1);
end

end


%% Algorithm functions

% Algorithm including optimization over friction exponent
function [X, Wfull, alphanew] = AlgAlpha(Wb,T,tol,n,na,na_idx,nu_idx,m,p_ip,b,b_full,d,MS1,MS2,data,dataFull,J,lowerB1,upperB1,lowerB2,upperB2,X0_1,regroup,Ginv,includeOffset,includeConstraints,constraintVariant,options)

% Initialize friction exponents
alphanew = ones(n,1);
alphaold = zeros(n,1);

    % start while loop for friction exponents alpha
    while(norm(alphaold-alphanew)>tol)
        
        % compute observation matrix for motor and friction dynamics 
        [Wfm,~,Wf] = Wextension(alphanew,n,na,na_idx,nu_idx,J,m,data,dataFull,includeOffset,includeConstraints,constraintVariant);

        % put together full observation matrix
        Wfull = [Wb Wfm];

        % initialize least squares function
        fun1 = @(X1) optimFun(X1, W, T);

        % formulate least squares problem for all parameters
        problem1 = createOptimProblem('fmincon','x0',X0_1,'objective',fun1,'lb',lowerB1,'ub',upperB1,'nonlcon',@(X1) lmiconDet(X1,n,Ginv,b,b_full,d,regroup),'options',options);
                
        % choose multi start option to increase the search space
        gs = MultiStart;

        % start first optimization
        X = run(gs,problem1,MS1);
                
        % set initial condition to identified inertial parameters
        X0_1 = X;       
    
        % check if regrouping is wanted 
        if regroup
    
        % case with regrouping
    
            % check if offset is a parameter
            if inlcludeOffset
                % compute estimated friction force/torque including offset
                Fest = (T - W*[X(1:b+n)' zeros(1,3*n)]');
                % get intial conditions for friction and friction exponent optimization problem
                X0_2 = [alphanew; X(b+n+1:end)];
            else
                % compute estimated friction force/torque without offset 
                Fest = (T - W*[X(1:b+n)' zeros(1,2*n)]');
                % get intial conditions for friction and friction exponent optimization problem
                X0_2 = [alphanew; X(b+n+1:end)];
            end
    
        else
    
         % case without regrouping
    
            % check if offset is a parameter
            if includeOffset
                % compute estimated friction force/torque including offset
                Fest = (T - W*[X(1:p_ip+n)' zeros(1,3*n)]');
                % get intial conditions for friction and friction exponent optimization problem
                X0_2 = [alphanew; X(11*n+1:end)];
            else
                % compute estimated friction force/torque without offset
                Fest = (T - W*[X(1:p_ip+n)' zeros(1,2*n)]');
                % get intial conditions for friction and friction exponent optimization problem
                X0_2 = [alphanew; X(11*n+1:end)];
            end
    
        end
        
        % update alpha for while loop check
        alphaold = alphanew;
    
        fun2 = @(X2) L2fric(Wf,X2,n,na,na_idx,nu_idx,m,data,dataFull,J,Fest,includeOffset,includeConstraints,constraintVariant);
        
        % formulate least squares problem for friction exponents and parameters
        problem2 = createOptimProblem('fmincon','x0',X0_2,'objective',fun2,'lb',lowerB2,'ub',upperB2,'options',options);
       
        % choose multi start option to increase the search space
        gs = MultiStart;
    
        % start second optimization
        X2 = run(gs,problem2,MS2);
    
        % update new friction exponents
        alphanew = X2(1:n);
    
        % check if offset is a parameter
        if includeOffset
            % update friction and offset parameters
            X(end-3*n+1:end) = X2(n+1:end);
        else
            % update friction parameters
            X(end-2*n+1:end) = X2(n+1:end);
        end
    
    end

end


% ------------------------------------------------------------------------


% Algorithm without optimimzation over friction exponent
function [X, Wfull] = AlgNoAlpha(Wb,T,n,na,na_idx,nu_idx,m,b,b_full,d,MS1,data,dataFull,J,lowerB1,upperB1,X0_1,regroup,Ginv,includeOffset,includeConstraints,constraintVariant,options)
        
     % compute observation matrix for motor and friction dynamics 
     [Wfm,~,~] = Wextension(ones(n,1),n,na,na_idx,nu_idx,J,m,data,dataFull,includeOffset,includeConstraints,constraintVariant);

     % put together full observation matrix
     Wfull = [Wb Wfm];

     % initialize least squares function
     fun1 = @(X1) optimFun(X1, W, T);

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

% ------------------------------------------------------------------------

% least squares function for friction parameters and exponent estimation
function [fmin, gradf] = L2fric(Wf,X,nLinks,na,na_idx,nu_idx,mData,data,dataFull,J,Fest,includeOffset,includeConstraints,constraintVariant)

    % get friction exponents
    alpha = X(1:nLinks);

    % get friction and offset parameters
    phiF = X(nLinks+1:end);

    % least squares function
    fmin = norm(Wf*phiF - Fest)^2;
    
    % gradient of least squares function
    if nargout > 1
        gradf = [2*L2FricGradAlpha(alpha,nLinks,na,na_idx,nu_idx,mData,data,dataFull,J,phiF,includeOffset,includeConstraints,constraintVariant).'*(Wf*phiF - Fest);
                 2*Wf.'*(Wf*phiF - Fest)];
    end

end

% ------------------------------------------------------------------------

% gradient of friction expnonent part of friction observation matrix
function fricGradAlpha = L2FricGradAlpha(alpha,nLinks,na,na_idx,nu_idx,mData,data,dataFull,J,phiF,includeOffset,includeConstraints,constraintVariant)

    % preallocate gradient matrix
    fricGradAlpha = zeros(nLinks*mData,nLinks);
    fricGrad_ = zeros(nLinks,nLinks);

    % small value to avoid log(0)
    log_eps = 0.000001;

    % compute gradients over all data point
    for i = 1:mData

        % compute gradient over each joint
        for j = 1:nLinks

            % check if offset is included
            if includeOffset
                fricGrad_(j,j) = [0 log(log_eps+abs(data(j,i + mData)))*sign(data(j,i + mData))*abs(data(j,i + mData))^alpha(j) 0]*phiF(3*(j-1)+1:3*j);
            else
                fricGrad_(j,j) = [0 sign(data(j,i + mData))*log(log_eps+abs(data(j,i + mData)))*abs(data(j,i + mData))^alpha(j)]*phiF(2*(j-1)+1:2*j);
            end

        end

        % assign gradient of one datapoint and all joints in the gradient matrix
        fricGradAlpha((i-1)*nLinks+1:i*nLinks,:) = fricGrad_;

    end

    % check if constraints are included
    if includeConstraints
        
        if constraintVariant == 1
        % general case for all kinds of constraints
        
            for i = 1:mData
        
                % compute constraint projection matrix
                J_temp = J(dataFull(:,i));
                K = eye(nLinks) - J_temp(:,sort([na_idx nu_idx])).'*pinv(J_temp(:,sort([na_idx nu_idx])).');
        
                % compute gradient of constrained observation matrix for friction exponent part
                fricGradAlpha((i-1)*nLinks+1:i*nLinks,:) = K*fricGradAlpha((i-1)*nLinks+1:i*nLinks,:);
            end
        
        elseif constraintVariant == 2
        % special case for fully actuated systems
        
            fricGradAlpha_temp = [];
        
            for i = 1:mData
        
                % compute constraint projection matrix
                J_temp = J(dataFull(:,i));
                P = [eye(na); -J_temp(:,nu_idx)\J_temp(:,na_idx)];
        
                % compute gradient of constrained observation matrix for friction exponent part on actuated joint level
                fricGradAlpha_ = P.'*fricGradAlpha((i-1)*nLinks+1:i*nLinks,:);
                fricGradAlpha_temp = [fricGradAlpha_temp;fricGradAlpha_];
            end
        
            % assign temporary gradient variable to gradient
            fricGradAlpha = fricGradAlpha_temp;
        
        end
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